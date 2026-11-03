"""Read-level amplicon QC reports, written when --qc-prefix is set."""

import csv
import re
from collections import Counter, defaultdict
from typing import NamedTuple, Optional

import numpy as np
import pysam
from primalbedtools.amplicons import Amplicon
from primalbedtools.bedfiles import BedLine

from align_trim.console import console

# reasons a fragment (single-end read or read pair) is removed, in pipeline order; reported as removed_<reason>
PRE_ASSIGNMENT_REASONS = [
    "unmapped",
    "not_proper_pair",
    "missing_mate",
    "low_mapq",
    "no_alignment",
]
AMPLICON_REASONS = [
    "no_primer_match",
    "partial_in_overlap",
    "mismatched_primers",
    "trim_failed",
    "primer_only",
    "not_full_length",
    "by_normalise",
]
REASONS = PRE_ASSIGNMENT_REASONS + AMPLICON_REASONS

PAIR_TYPES = [
    "expected",
    "overlap_product",
    "partial_in_overlap",
    "long_same_pool",
    "long_cross_pool",
    "cross_chrom",
    "other",
    "two_left_primers",
    "two_right_primers",
]

# split read join types, most to least notable
SPLIT_TYPES = ["amplicon_chimera", "within_amplicon", "low_mapq"]

# default --primer-match-threshold: max distance (bp) between a read end and a primer 5' end (or amplicon end) for
# the end to count as at that primer, used for primer lookup and every QC end check
DEFAULT_PRIMER_MATCH_THRESHOLD = 5

# bases of soft-clipped barcode/adapter allowed for, beyond the primer itself, when scanning read ends for primer k-mers
SCAN_ADAPTER_FLANK = 35

# overlap concordance thresholds
OVERLAP_MIN_DEPTH = 10
OVERLAP_AF_THRESHOLD = 0.2

# split alignment breakpoint cluster bin size (bp)
BREAKPOINT_BIN = 10

# length of the primer 3' k-mer searched for at unexplained read ends
SAME_DIRECTION_K = 12

# CIGAR op codes which are matches (M, =, X)
MATCH_OPS = (0, 7, 8)
CIGAR_OPS = "MIDNSHP=XB"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# overlap count columns
COUNT_COLS = ["A", "C", "G", "T", "N", "DEL", "INS"]
BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
N_IDX, DEL_IDX, INS_IDX = 4, 5, 6
ALLELE_COLS = [0, 1, 2, 3, DEL_IDX]  # alleles contributing to depth

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def classify_primer_pair(
    p1: BedLine,
    p2: BedLine,
    start=None,
    end=None,
    tolerance=DEFAULT_PRIMER_MATCH_THRESHOLD,
) -> str:
    """Classify a (left, right) primer pair by what it means in a tiling scheme.

    start/end are the fragment's reference span. When given, an overlap product must start and end within
    `tolerance` of the primer 5' ends, otherwise it is a partial piece of a normal amplicon which happens to
    lie in the overlap.
    """
    if p1.chrom != p2.chrom:
        return "cross_chrom"
    k = p2.amplicon_number - p1.amplicon_number
    if k == 0:
        return "expected"
    if k == -1:
        # LEFT_{n+1} / RIGHT_n, short product from the overlap region
        if start is not None and (
            abs(start - p1.start) > tolerance or abs(end - p2.end) > tolerance
        ):
            return "partial_in_overlap"
        return "overlap_product"
    if k >= 1:
        # LEFT_n / RIGHT_{n+k}, the intervening amplicon(s) were skipped
        return "long_same_pool" if p1.pool == p2.pool else "long_cross_pool"
    return "other"


class Part(NamedTuple):
    """One aligned part of a read. Query coordinates are in read orientation."""

    chrom: str
    ref_start: int
    ref_end: int
    strand: str
    mapq: int
    q_start: int
    q_end: int


def make_part(chrom, ref_start, strand, mapq, ops) -> Part:
    """Build a Part from a list of (length, op) CIGAR tuples where op is a CIGAR character."""
    ref_len = sum(length for length, op in ops if op in "MDN=X")
    q_len = sum(length for length, op in ops if op in "MIS=XH")
    lead = 0
    for length, op in ops:
        if op not in "SH":
            break
        lead += length
    trail = 0
    for length, op in reversed(ops):
        if op not in "SH":
            break
        trail += length
    q_start, q_end = lead, q_len - trail
    if strand == "-":
        q_start, q_end = q_len - q_end, q_len - q_start
    return Part(chrom, ref_start, ref_start + ref_len, strand, mapq, q_start, q_end)


def parse_sa_tag(sa: str) -> list[Part]:
    """Parse a SAM SA tag (rname,pos,strand,CIGAR,mapQ,NM;...) into Parts."""
    parts = []
    for entry in sa.split(";"):
        if not entry:
            continue
        chrom, pos, strand, cigar, mapq, _ = entry.split(",")
        ops = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
        parts.append(make_part(chrom, int(pos) - 1, strand, int(mapq), ops))
    return parts


def primary_part(segment: pysam.AlignedSegment) -> Part:
    ops = [(length, CIGAR_OPS[op]) for op, length in segment.cigartuples]  # type: ignore
    return make_part(
        segment.reference_name,
        segment.reference_start,
        "-" if segment.is_reverse else "+",
        segment.mapping_quality,
        ops,
    )


def junctions(parts: list[Part]) -> list[tuple[Part, Part]]:
    """Consecutive pairs of parts in read order, each pair being one join in the read."""
    ordered = sorted(parts, key=lambda p: p.q_start)
    return list(zip(ordered, ordered[1:]))


def breakpoints(a: Part, b: Part) -> tuple[int, int]:
    """0-based interbase reference coordinates of the break at the end of part a and the start of part b."""
    return (
        a.ref_end if a.strand == "+" else a.ref_start,
        b.ref_start if b.strand == "+" else b.ref_end,
    )


def breakpoint_positions(a: Part, b: Part) -> tuple[int, int]:
    """1-based positions of the last aligned base before the break (in part a) and the first after it (in part b)."""
    return (
        a.ref_end if a.strand == "+" else a.ref_start + 1,
        b.ref_start + 1 if b.strand == "+" else b.ref_end,
    )


def join_orientation(a: Part, b: Part) -> str:
    """inverted: parts on opposite strands; collinear: same strand and in reference order; else non_collinear."""
    if a.strand != b.strand:
        return "inverted"
    if a.chrom == b.chrom and (
        b.ref_start >= a.ref_start if a.strand == "+" else b.ref_start <= a.ref_start
    ):
        return "collinear"
    return "non_collinear"


def count_bases_in_window(segment, start, end, counts, seen=None):
    """Add a segment's bases over reference [start, end) to counts.

    counts has one row per window position and columns A C G T N DEL INS. Insertions are
    attributed to the reference position before them. Positions flagged in `seen` (already
    counted for the template's other mate) are skipped.

    Returns a boolean array of the window positions this segment counted.
    """
    covered = np.zeros(end - start, dtype=bool)
    seq = segment.query_sequence
    if (
        seq is None
        or segment.reference_end is None
        or segment.reference_end <= start
        or segment.reference_start >= end
    ):
        return covered

    def skip(i):
        return seen is not None and seen[i]

    rpos = segment.reference_start
    qpos = 0
    for op, length in segment.cigartuples:
        if op in MATCH_OPS:
            for p in range(max(rpos, start), min(rpos + length, end)):
                i = p - start
                if skip(i):
                    continue
                counts[i, BASE_INDEX.get(seq[qpos + p - rpos], N_IDX)] += 1
                covered[i] = True
            rpos += length
            qpos += length
        elif op == 1:
            i = rpos - 1 - start
            if 0 <= i < end - start and not skip(i):
                counts[i, INS_IDX] += 1
            qpos += length
        elif op in (2, 3):
            for p in range(max(rpos, start), min(rpos + length, end)):
                i = p - start
                if skip(i):
                    continue
                counts[i, DEL_IDX] += 1
                covered[i] = True
            rpos += length
        elif op == 4:
            qpos += length
        if rpos > end:
            break
    return covered


def overlap_stats(left, right):
    """Per-position depth, AF difference, evaluability and major alleles for one overlap."""
    left_depth = left[:, ALLELE_COLS].sum(axis=1)
    right_depth = right[:, ALLELE_COLS].sum(axis=1)
    cols = ALLELE_COLS + [INS_IDX]
    with np.errstate(divide="ignore", invalid="ignore"):
        left_af = np.nan_to_num(left[:, cols] / left_depth[:, None])
        right_af = np.nan_to_num(right[:, cols] / right_depth[:, None])
    af_diff = np.abs(left_af - right_af).max(axis=1)
    evaluable = (left_depth >= OVERLAP_MIN_DEPTH) & (right_depth >= OVERLAP_MIN_DEPTH)
    left_major = np.array([COUNT_COLS[c] for c in ALLELE_COLS])[
        left[:, ALLELE_COLS].argmax(axis=1)
    ]
    right_major = np.array([COUNT_COLS[c] for c in ALLELE_COLS])[
        right[:, ALLELE_COLS].argmax(axis=1)
    ]
    return left_depth, right_depth, af_diff, evaluable, left_major, right_major


def base_frequencies(counts):
    """Relative frequency of each base (A, C, G, T, DEL) at one position, most common first, as `C=0.95,T=0.05`.

    Bases below 0.01 are left out.
    """
    depth = counts[ALLELE_COLS].sum()
    freqs = sorted(
        ((COUNT_COLS[c], counts[c] / depth) for c in ALLELE_COLS),
        key=lambda x: -x[1],
    )
    return ",".join(f"{base}={freq:.2f}" for base, freq in freqs if freq >= 0.005)


def _write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class QC:
    """Accumulates the --qc-prefix reports during the single pass over the input."""

    def __init__(
        self,
        amplicon_list: list[Amplicon],
        lookup: dict,
        bedlines: list[BedLine],
        min_mapq: int,
        padding: int,
    ):
        # avoid a circular import, main imports this module
        from align_trim.main import find_primer_with_lookup

        self._find_primer = find_primer_with_lookup
        self.lookup = lookup
        self.min_mapq = min_mapq
        self.padding = padding
        self.amplicons = {(a.chrom, a.amplicon_number): a for a in amplicon_list}

        # kept/removed fragment counts per amplicon, pre-seeded so amplicons with no fragments still get a row,
        # plus fragments removed before they could be assigned to an amplicon
        self.fates = defaultdict(Counter)
        for key in sorted(self.amplicons):
            self.fates[key]
        self.unassigned = Counter()

        # primer pairs, keyed by (left BedLine, right BedLine, type)
        self.pairs = Counter()

        # split reads: fragments checked, fragments by most notable join type, and joins
        self.split_checked = 0
        self.split_fragments = Counter()
        self.joins = Counter()
        outer_ends = defaultdict(set)
        for a in amplicon_list:
            outer_ends[a.chrom].update((a.amplicon_start, a.amplicon_end))
        self.outer_ends = {c: np.array(sorted(e)) for c, e in outer_ends.items()}

        # overlap windows between adjacent amplicons: [cov_start(n+1), cov_end(n))
        self.windows = []
        self.windows_by_amp = defaultdict(list)
        by_chrom = defaultdict(dict)
        for a in amplicon_list:
            by_chrom[a.chrom][a.amplicon_number] = a
        for chrom, amps in sorted(by_chrom.items()):
            for n in sorted(amps):
                if n + 1 not in amps:
                    continue
                start, end = amps[n + 1].coverage_start, amps[n].coverage_end
                if end <= start:
                    continue
                wid = len(self.windows)
                self.windows.append(
                    {
                        "chrom": chrom,
                        "amplicon_left": n,
                        "amplicon_right": n + 1,
                        "start": start,
                        "end": end,
                        "left": np.zeros((end - start, len(COUNT_COLS)), dtype=int),
                        "right": np.zeros((end - start, len(COUNT_COLS)), dtype=int),
                    }
                )
                self.windows_by_amp[(chrom, n)].append((wid, "left"))
                self.windows_by_amp[(chrom, n + 1)].append((wid, "right"))

        # primer 3' k-mers as they appear on the reference-oriented read at each end
        self.left_end_kmers = {}
        self.right_end_kmers = {}
        for bl in bedlines:
            if not bl.sequence or len(bl.sequence) < SAME_DIRECTION_K:
                continue
            kmer = bl.sequence[-SAME_DIRECTION_K:].upper()
            self.left_end_kmers.setdefault(kmer, bl)
            self.right_end_kmers.setdefault(revcomp(kmer), bl)
        self.scan_window = max((len(bl.sequence or "") for bl in bedlines), default=0)
        self.scan_window += SCAN_ADAPTER_FLANK

    @classmethod
    def from_args(cls, args, amplicon_list, lookup, bedlines) -> Optional["QC"]:
        if not getattr(args, "qc_prefix", None):
            return None
        return cls(
            amplicon_list,
            lookup,
            bedlines,
            min_mapq=args.min_mapq,
            padding=args.primer_match_threshold,
        )

    # fragment fate: "kept" or one of REASONS
    def fate(self, fate, chrom=None, amplicon=None):
        if amplicon is None:
            self.unassigned[fate] += 1
        else:
            self.fates[(chrom, amplicon)][fate] += 1

    def kept(self, segments, chrom, amplicon, correctly_paired):
        self.fate("kept", chrom, amplicon)
        if correctly_paired:
            self._count_overlaps(segments, chrom, amplicon)

    # primer pairs
    def pair(self, p1, p2, pair_type=None, start=None, end=None):
        """Count a primer pair and return its type."""
        pair_type = pair_type or classify_primer_pair(p1, p2, start, end, self.padding)
        self.pairs[(p1, p2, pair_type)] += 1
        return pair_type

    def same_direction(self, left_seg, right_seg, p1, p2):
        """Look for a same-direction primer at read ends not explained by the assigned primers.

        Records a two_left_primers or two_right_primers pair when found. Never changes what happens to the read.
        """
        left_ok = p1 and abs(left_seg.reference_start - p1.start) <= self.padding
        right_ok = p2 and abs(right_seg.reference_end - p2.end) <= self.padding
        if left_ok and right_ok:
            return
        left_hit = None if left_ok else self._scan_end(left_seg, left_end=True)
        right_hit = None if right_ok else self._scan_end(right_seg, left_end=False)
        left = p1 if left_ok else left_hit
        right = p2 if right_ok else right_hit
        if not left or not right:
            return
        if left.strand == "+" and right.strand == "+":
            self.pair(left, right, "two_left_primers")
        elif left.strand == "-" and right.strand == "-":
            self.pair(left, right, "two_right_primers")

    def _scan_end(self, segment, left_end):
        seq = segment.query_sequence
        if not seq:
            return None
        k = SAME_DIRECTION_K
        if left_end:
            region = seq[: self.scan_window]
            for i in range(len(region) - k + 1):
                hit = self.left_end_kmers.get(region[i : i + k])
                if hit:
                    return hit
        else:
            region = seq[-self.scan_window :]
            for i in range(len(region) - k, -1, -1):
                hit = self.right_end_kmers.get(region[i : i + k])
                if hit:
                    return hit
        return None

    # split alignments
    def _near_amplicon_end(self, chrom, pos):
        ends = self.outer_ends.get(chrom)
        if ends is None:
            return False
        i = np.searchsorted(ends, pos)
        return any(
            abs(int(ends[j]) - pos) <= self.padding
            for j in (i - 1, i)
            if 0 <= j < len(ends)
        )

    def _part_amplicon(self, part):
        if part.chrom not in self.lookup:
            return "NA"
        if not 0 <= part.ref_start < self.lookup[part.chrom].shape[1]:
            return "NA"
        p = self._find_primer(self.lookup, part.ref_start, "+", part.chrom)
        return p.amplicon_number if p else "NA"

    def split(self, segments):
        """Count the joins in a fragment's split alignments (SA tags of each of its reads).

        A join seen in both reads of a pair is counted once for the fragment.
        """
        self.split_checked += 1
        rows = set()
        for segment in segments:
            if not segment.has_tag("SA"):
                continue
            parts = [primary_part(segment)] + parse_sa_tag(segment.get_tag("SA"))  # type: ignore
            for a, b in junctions(parts):
                if a.mapq < self.min_mapq or b.mapq < self.min_mapq:
                    join_type = "low_mapq"
                elif all(
                    self._near_amplicon_end(part.chrom, pos)
                    for part, pos in zip((a, b), breakpoints(a, b))
                ):
                    join_type = "amplicon_chimera"
                else:
                    join_type = "within_amplicon"
                # each side is (chrom, 10 bp window start, amplicon), lower coordinate side first
                sides = sorted(
                    (
                        (
                            part.chrom,
                            (pos - 1) // BREAKPOINT_BIN * BREAKPOINT_BIN + 1,
                            self._part_amplicon(part),
                        )
                        for part, pos in zip((a, b), breakpoint_positions(a, b))
                    ),
                    key=lambda side: side[:2],
                )
                rows.add((*sides[0], *sides[1], join_type, join_orientation(a, b)))
        if not rows:
            return
        self.split_fragments[min((r[6] for r in rows), key=SPLIT_TYPES.index)] += 1
        for row in rows:
            self.joins[row] += 1

    # overlap concordance
    def _count_overlaps(self, segments, chrom, amplicon):
        for wid, side in self.windows_by_amp.get((chrom, amplicon), []):
            w = self.windows[wid]
            seen = None
            for seg in segments:
                covered = count_bases_in_window(
                    seg, w["start"], w["end"], w[side], seen
                )
                seen = covered if seen is None else seen | covered

    # output
    def write(self, prefix):
        console.heading("QC summary")
        self._write_summary(f"{prefix}.summary.tsv")
        self._write_amplicons(f"{prefix}.amplicons.tsv")
        self._write_primer_pairs(f"{prefix}.primer_pairs.tsv")
        self._write_split_reads(f"{prefix}.split_reads.tsv")
        self._write_overlaps(f"{prefix}.overlaps.tsv")

    def _write_summary(self, path):
        totals = Counter(self.unassigned)
        for counts in self.fates.values():
            totals.update(counts)
        fragments = sum(totals.values())
        row = {"fragments": fragments, "kept": totals["kept"]}
        row.update({f"removed_{r}": totals[r] for r in REASONS})
        _write_tsv(path, list(row), [row])

        removed = Counter({r: totals[r] for r in REASONS})
        most_removed = ", ".join(
            f"{r.replace('_', ' ')} {n} ({_pct(n, fragments)})"
            for r, n in removed.most_common(3)
            if n
        )
        console.prefixed(
            "QC:",
            f"{fragments} fragments, {totals['kept']} kept ({_pct(totals['kept'], fragments)})"
            + (f"; most removed: {most_removed}" if most_removed else ""),
        )

    def _write_amplicons(self, path):
        fields = ["chrom", "amplicon", "pool", "start", "end", "fragments", "kept"]
        fields += [f"removed_{r}" for r in AMPLICON_REASONS]
        rows = []
        for (chrom, number), counts in self.fates.items():
            amp = self.amplicons[(chrom, number)]
            row = {
                "chrom": chrom,
                "amplicon": number,
                "pool": amp.pool,
                "start": amp.amplicon_start + 1,
                "end": amp.amplicon_end,
                "fragments": sum(counts.values()),
                "kept": counts["kept"],
            }
            row.update({f"removed_{r}": counts[r] for r in AMPLICON_REASONS})
            rows.append(row)
        _write_tsv(path, fields, rows)

        no_fragments = sum(1 for r in rows if r["fragments"] == 0)
        none_kept = sum(1 for r in rows if r["fragments"] and not r["kept"])
        console.prefixed(
            "QC:",
            f"amplicons: {len(rows)} in scheme, {no_fragments} with no fragments, "
            f"{none_kept} with fragments but none kept",
        )

    def _write_primer_pairs(self, path):
        fields = [
            "left_chrom",
            "left_primer",
            "left_amplicon",
            "right_chrom",
            "right_primer",
            "right_amplicon",
            "type",
            "fragments",
        ]
        rows = [
            {
                "left_chrom": p1.chrom,
                "left_primer": p1.primername,
                "left_amplicon": p1.amplicon_number,
                "right_chrom": p2.chrom,
                "right_primer": p2.primername,
                "right_amplicon": p2.amplicon_number,
                "type": pair_type,
                "fragments": n,
            }
            for (p1, p2, pair_type), n in self.pairs.items()
        ]
        rows.sort(
            key=lambda r: (
                r["left_chrom"],
                r["left_amplicon"],
                r["right_chrom"],
                r["right_amplicon"],
                PAIR_TYPES.index(r["type"]),
            )
        )
        _write_tsv(path, fields, rows)

        totals = Counter()
        for (_, _, pair_type), n in self.pairs.items():
            totals[pair_type] += n

        t = totals
        console.prefixed(
            "QC:",
            f"primer pairs: same amplicon (expected) {t['expected']}; "
            f"short products of an overlap region {t['overlap_product']}; "
            f"fragments inside an overlap with ends not at its primers {t['partial_in_overlap']}; "
            f"long products skipping a primer site: same pool {t['long_same_pool']}, "
            f"cross pool {t['long_cross_pool']}; "
            f"between chromosomes {t['cross_chrom']}; other {t['other']}; "
            f"also misprimed with two left primers {t['two_left_primers']}, "
            f"two right primers {t['two_right_primers']}",
        )

    def _write_split_reads(self, path):
        fields = [
            "chrom",
            "pos",
            "amplicon",
            "joined_chrom",
            "joined_pos",
            "joined_amplicon",
            "type",
            "orientation",
            "fragments",
        ]
        rows = [dict(zip(fields, (*key, n))) for key, n in self.joins.items()]
        rows.sort(key=lambda r: -r["fragments"])
        _write_tsv(path, fields, rows)

        split = sum(self.split_fragments.values())
        chimeras = self.split_fragments["amplicon_chimera"]
        console.prefixed(
            "QC:",
            f"split reads: {split} of {self.split_checked} fragments ({_pct(split, self.split_checked)}); "
            f"amplicon chimeras {chimeras} ({_pct(chimeras, self.split_checked)}), "
            f"within amplicon {self.split_fragments['within_amplicon']}, "
            f"low mapq {self.split_fragments['low_mapq']}",
        )

    def _write_overlaps(self, path):
        fields = [
            "chrom",
            "amplicon",
            "next_amplicon",
            "start",
            "end",
            "length",
            "compared_positions",
            "amplicon_median_depth",
            "next_amplicon_median_depth",
            "max_freq_difference",
            "max_freq_difference_pos",
            "differing_positions",
            "majority_base_differences",
            "majority_base_difference_details",
        ]
        rows = []
        for w in self.windows:
            depth, next_depth, freq_diff, compared, major, next_major = overlap_stats(
                w["left"], w["right"]
            )
            row = {
                "chrom": w["chrom"],
                "amplicon": w["amplicon_left"],
                "next_amplicon": w["amplicon_right"],
                "start": w["start"] + 1,
                "end": w["end"],
                "length": w["end"] - w["start"],
                "compared_positions": int(compared.sum()),
                "amplicon_median_depth": round(float(np.median(depth)), 2),
                "next_amplicon_median_depth": round(float(np.median(next_depth)), 2),
            }
            if compared.any():
                i = int(np.where(compared, freq_diff, -1).argmax())
                differs = compared & (major != next_major)
                row["max_freq_difference"] = round(float(freq_diff[i]), 4)
                row["max_freq_difference_pos"] = w["start"] + i + 1
                row["differing_positions"] = int(
                    (compared & (freq_diff >= OVERLAP_AF_THRESHOLD)).sum()
                )
                row["majority_base_differences"] = int(differs.sum())
                row["majority_base_difference_details"] = ";".join(
                    f"{w['start'] + i + 1}:{base_frequencies(w['left'][i])}/{base_frequencies(w['right'][i])}"
                    for i in np.flatnonzero(differs)
                )
            else:
                for field in fields[9:]:
                    row[field] = "NA"
            rows.append(row)
        _write_tsv(path, fields, rows)

        compared_rows = [r for r in rows if r["compared_positions"]]
        differing = sum(1 for r in compared_rows if r["differing_positions"])
        mismatched = [r for r in compared_rows if r["majority_base_differences"]]
        console.prefixed(
            "QC:",
            f"amplicon overlaps: {len(compared_rows)} of {len(rows)} had enough depth to compare; "
            f"base frequencies differ in {differing}, majority base differs in {len(mismatched)}",
        )
        if mismatched:
            listed = ", ".join(
                f"{r['amplicon']}/{r['next_amplicon']} ({r['majority_base_differences']} "
                f"position{'s' if r['majority_base_differences'] > 1 else ''})"
                for r in mismatched
            )
            warning = (
                f"majority base differs between overlapping amplicons {listed}, "
                f"which can indicate contamination; see {path}"
            )
            console.defer_warning(warning)
            return warning
        return None


def _pct(n, total):
    return f"{(n / total) * 100:.2f}%" if total else "0.00%"

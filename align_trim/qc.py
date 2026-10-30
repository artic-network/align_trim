"""Read-level amplicon QC reports, written when --qc-prefix is set."""

import csv
import re
import sys
from collections import Counter, defaultdict
from typing import NamedTuple, Optional

import numpy as np
import pysam
from primalbedtools.amplicons import Amplicon
from primalbedtools.bedfiles import BedLine

# read fates, in report column order
FATES = [
    "kept",
    "normalised_out",
    "not_full_length",
    "fully_masked",
    "softmask_error",
    "incorrectly_paired",
    "no_primer",
    "low_mapq",
    "no_mapping_data",
    "unmapped",
    "supplementary",
    "secondary",
    "missing_mate",
]
UNASSIGNED = "unassigned"

PAIR_CATEGORIES = [
    "correct",
    "short_overlap",
    "long_same_pool",
    "long_cross_pool",
    "cross_chrom",
    "other",
    "left_left",
    "right_right",
]

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


def classify_primer_pair(p1: BedLine, p2: BedLine) -> str:
    """Classify a (left, right) primer pair by what it means in a tiling scheme."""
    if p1.chrom != p2.chrom:
        return "cross_chrom"
    k = p2.amplicon_number - p1.amplicon_number
    if k == 0:
        return "correct"
    if k == -1:
        # LEFT_{n+1} / RIGHT_n, short product from the overlap region
        return "short_overlap"
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


def junctions(parts: list[Part]) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    """Reference breakpoints between consecutive parts in read order."""
    ordered = sorted(parts, key=lambda p: p.q_start)
    out = []
    for a, b in zip(ordered, ordered[1:]):
        bp_a = a.ref_end if a.strand == "+" else a.ref_start
        bp_b = b.ref_start if b.strand == "+" else b.ref_end
        out.append(((a.chrom, bp_a), (b.chrom, bp_b)))
    return out


def split_orientation(primary: Part, parts: list[Part]) -> str:
    if any(p.strand != primary.strand for p in parts):
        return "inverted"
    ordered = sorted(parts, key=lambda p: p.q_start)
    if len({p.chrom for p in ordered}) == 1:
        starts = [p.ref_start for p in ordered]
        if primary.strand == "-":
            starts = starts[::-1]
        if starts == sorted(starts):
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

        # read fate, pre-seeded so amplicons with no reads still get a row
        self.fates = defaultdict(Counter)
        for key in sorted(self.amplicons):
            self.fates[key]

        # primer pairs, keyed by (left BedLine, right BedLine, category)
        self.pairs = Counter()

        # split alignments
        self.split_reads = 0
        self.split_classes = Counter()
        self.splits = Counter()
        self.breakpoints = Counter()
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
        self.scan_window += padding

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

    # read fate
    def fate(self, fate, chrom=None, amplicon=None):
        if amplicon is None:
            key = (chrom or "*", UNASSIGNED)
        else:
            key = (chrom, amplicon)
        self.fates[key][fate] += 1

    def kept(self, segments, chrom, amplicon, correctly_paired):
        self.fate("kept", chrom, amplicon)
        if correctly_paired:
            self._count_overlaps(segments, chrom, amplicon)

    # primer pairs
    def pair(self, p1, p2, category=None):
        self.pairs[(p1, p2, category or classify_primer_pair(p1, p2))] += 1

    def same_direction(self, left_seg, right_seg, p1, p2):
        """Look for a same-direction primer at read ends not explained by the assigned primers.

        Records a left_left or right_right pair when found. Never changes what happens to the read.
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
            self.pair(left, right, "left_left")
        elif left.strand == "-" and right.strand == "-":
            self.pair(left, right, "right_right")

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

    def split(self, segment):
        self.split_reads += 1
        if not segment.has_tag("SA"):
            return
        primary = primary_part(segment)
        parts = [primary] + parse_sa_tag(segment.get_tag("SA"))  # type: ignore
        joins = junctions(parts)
        if any(p.mapq < self.min_mapq for p in parts):
            cls = "low_confidence"
        elif all(
            self._near_amplicon_end(chrom, pos) for join in joins for chrom, pos in join
        ):
            cls = "amplicon_junction"
        else:
            cls = "internal_breakpoint"
        orientation = split_orientation(primary, parts)
        self.split_classes[cls] += 1

        def amp_key(x):
            return (x[0], x[1] if isinstance(x[1], int) else float("inf"))

        ordered = sorted(parts, key=lambda p: p.q_start)
        a, b = sorted(
            [
                (ordered[0].chrom, self._part_amplicon(ordered[0])),
                (ordered[1].chrom, self._part_amplicon(ordered[1])),
            ],
            key=amp_key,
        )
        self.splits[(cls, orientation, *a, *b)] += 1

        clusters = set()
        for join in joins:
            bins = sorted(
                (chrom, pos // BREAKPOINT_BIN * BREAKPOINT_BIN) for chrom, pos in join
            )
            clusters.add((*bins[0], *bins[1], cls, orientation))
        for cluster in clusters:
            self.breakpoints[cluster] += 1

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
        self._write_fates(f"{prefix}.read_fate.tsv")
        self._write_pairs(f"{prefix}.primer_pairs.tsv")
        self._write_splits(
            f"{prefix}.split_alignments.tsv", f"{prefix}.split_breakpoints.tsv"
        )
        self._write_overlaps(
            f"{prefix}.overlap_positions.tsv", f"{prefix}.overlap_summary.tsv"
        )

    def _write_fates(self, path):
        assigned = [k for k in self.fates if k[1] != UNASSIGNED]
        unassigned = sorted(k for k in self.fates if k[1] == UNASSIGNED)
        rows = []
        for key in assigned + unassigned:
            counts = self.fates[key]
            amp = self.amplicons.get(key)
            row = {
                "chrom": key[0],
                "amplicon": key[1],
                "pool": amp.pool if amp else "NA",
                "total": sum(counts.values()),
            }
            row.update({f: counts[f] for f in FATES})
            rows.append(row)
        _write_tsv(path, ["chrom", "amplicon", "pool", "total"] + FATES, rows)

        total = sum(r["total"] for r in rows)
        kept = sum(r["kept"] for r in rows)
        normalised_out = sum(r["normalised_out"] for r in rows)
        amp_rows = [r for r in rows if r["amplicon"] != UNASSIGNED]
        no_reads = sum(1 for r in amp_rows if r["total"] == 0)
        none_kept = sum(1 for r in amp_rows if r["total"] > 0 and r["kept"] == 0)
        pct = (kept / total) * 100 if total else 0
        print(
            f"Read fate: {total} templates, {kept} kept ({pct:.2f}%), {normalised_out} normalised out; "
            f"{no_reads} amplicons with 0 assigned, {none_kept} with assigned but 0 kept",
            file=sys.stderr,
        )

    def _write_pairs(self, path):
        rows = []
        for (p1, p2, category), count in self.pairs.items():
            rows.append(
                {
                    "chrom_left": p1.chrom,
                    "left_primer": p1.primername,
                    "left_amplicon": p1.amplicon_number,
                    "left_pool": p1.pool,
                    "chrom_right": p2.chrom,
                    "right_primer": p2.primername,
                    "right_amplicon": p2.amplicon_number,
                    "right_pool": p2.pool,
                    "category": category,
                    "count": count,
                }
            )
        rows.sort(
            key=lambda r: (
                r["chrom_left"],
                r["left_amplicon"],
                r["chrom_right"],
                r["right_amplicon"],
                r["category"],
            )
        )
        _write_tsv(path, list(rows[0]) if rows else ["category", "count"], rows)

        totals = Counter()
        for (_, _, category), count in self.pairs.items():
            totals[category] += count
        summary = " ".join(f"{c}={totals[c]}" for c in PAIR_CATEGORIES)
        print(f"Primer pairs: {summary}", file=sys.stderr)

    def _write_splits(self, path, breakpoint_path):
        fields = [
            "class",
            "orientation",
            "chrom_a",
            "amplicon_a",
            "chrom_b",
            "amplicon_b",
            "reads",
        ]
        rows = [dict(zip(fields, (*key, n))) for key, n in self.splits.items()]
        rows.sort(key=lambda r: -r["reads"])
        _write_tsv(path, fields, rows)

        fields = [
            "chrom_a",
            "bin_a",
            "chrom_b",
            "bin_b",
            "class",
            "orientation",
            "reads",
        ]
        rows = []
        for (
            chrom_a,
            bin_a,
            chrom_b,
            bin_b,
            cls,
            orientation,
        ), n in self.breakpoints.items():
            rows.append(
                {
                    "chrom_a": chrom_a,
                    "bin_a": bin_a + 1,
                    "chrom_b": chrom_b,
                    "bin_b": bin_b + 1,
                    "class": cls,
                    "orientation": orientation,
                    "reads": n,
                }
            )
        rows.sort(key=lambda r: -r["reads"])
        _write_tsv(breakpoint_path, fields, rows)

        with_sa = sum(self.split_classes.values())
        chimeric = self.split_classes["amplicon_junction"]
        pct = (chimeric / self.split_reads) * 100 if self.split_reads else 0
        print(
            f"Split alignments: {with_sa} of {self.split_reads} reads have SA tags; "
            f"chimera rate (amplicon_junction) {chimeric} ({pct:.2f}%); "
            f"internal_breakpoint {self.split_classes['internal_breakpoint']}; "
            f"low_confidence {self.split_classes['low_confidence']}",
            file=sys.stderr,
        )

    def _write_overlaps(self, path, summary_path):
        pos_fields = ["chrom", "amplicon_left", "amplicon_right", "pos", "left_depth"]
        pos_fields += [f"left_{c}" for c in COUNT_COLS] + ["right_depth"]
        pos_fields += [f"right_{c}" for c in COUNT_COLS]
        pos_fields += ["left_major", "right_major", "af_diff", "evaluable"]
        summary_fields = [
            "chrom",
            "amplicon_left",
            "amplicon_right",
            "start",
            "end",
            "length",
            "evaluable_positions",
            "median_left_depth",
            "median_right_depth",
            "max_af_diff",
            "max_af_diff_pos",
            "discordant_positions",
            "major_allele_mismatches",
        ]
        pos_rows = []
        summary_rows = []
        for w in self.windows:
            left, right = w["left"], w["right"]
            left_depth, right_depth, af_diff, evaluable, left_major, right_major = (
                overlap_stats(left, right)
            )
            for i in range(w["end"] - w["start"]):
                row = {
                    "chrom": w["chrom"],
                    "amplicon_left": w["amplicon_left"],
                    "amplicon_right": w["amplicon_right"],
                    "pos": w["start"] + i + 1,
                    "left_depth": int(left_depth[i]),
                    "right_depth": int(right_depth[i]),
                    "left_major": left_major[i],
                    "right_major": right_major[i],
                    "af_diff": round(float(af_diff[i]), 4),
                    "evaluable": bool(evaluable[i]),
                }
                for j, c in enumerate(COUNT_COLS):
                    row[f"left_{c}"] = int(left[i, j])
                    row[f"right_{c}"] = int(right[i, j])
                pos_rows.append(row)

            discordant = evaluable & (af_diff >= OVERLAP_AF_THRESHOLD)
            mismatched = evaluable & (left_major != right_major)
            if evaluable.any():
                masked = np.where(evaluable, af_diff, -1)
                max_i = int(masked.argmax())
                max_af_diff = round(float(af_diff[max_i]), 4)
                max_af_diff_pos = w["start"] + max_i + 1
            else:
                max_af_diff = max_af_diff_pos = "NA"
            summary_rows.append(
                {
                    "chrom": w["chrom"],
                    "amplicon_left": w["amplicon_left"],
                    "amplicon_right": w["amplicon_right"],
                    "start": w["start"] + 1,
                    "end": w["end"],
                    "length": w["end"] - w["start"],
                    "evaluable_positions": int(evaluable.sum()),
                    "median_left_depth": float(np.median(left_depth)),
                    "median_right_depth": float(np.median(right_depth)),
                    "max_af_diff": max_af_diff,
                    "max_af_diff_pos": max_af_diff_pos,
                    "discordant_positions": int(discordant.sum()),
                    "major_allele_mismatches": int(mismatched.sum()),
                }
            )
        _write_tsv(path, pos_fields, pos_rows)
        _write_tsv(summary_path, summary_fields, summary_rows)

        evaluable_overlaps = sum(1 for r in summary_rows if r["evaluable_positions"])
        discordant_overlaps = [r for r in summary_rows if r["discordant_positions"]]
        mismatch_overlaps = [
            f"{r['amplicon_left']}/{r['amplicon_right']}"
            for r in summary_rows
            if r["major_allele_mismatches"]
        ]
        listed = ",".join(mismatch_overlaps[:10])
        if len(mismatch_overlaps) > 10:
            listed += ",..."
        print(
            f"Overlap concordance: {evaluable_overlaps}/{len(summary_rows)} overlaps evaluable; "
            f"{len(discordant_overlaps)} with >=1 discordant position; "
            f"{len(mismatch_overlaps)} with major allele mismatch"
            + (f" ({listed})" if listed else ""),
            file=sys.stderr,
        )

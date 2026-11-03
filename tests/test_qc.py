import contextlib
import copy
import csv
import io
import pathlib
import random
import statistics
import tempfile
import unittest
from collections import Counter, namedtuple

import numpy as np
import pysam
from primalbedtools.amplicons import create_amplicons
from primalbedtools.bedfiles import merge_primers
from primalbedtools.scheme import Scheme

from align_trim.console import console
from align_trim.main import (
    create_primer_lookup,
    find_primer_with_lookup,
    fragment_spans,
    go,
    handle_segments,
    read_pair_generator,
    trim,
)
from align_trim.qc import (
    AMPLICON_REASONS,
    QC,
    base_frequencies,
    REASONS,
    breakpoint_positions,
    classify_primer_pair,
    count_bases_in_window,
    join_orientation,
    junctions,
    overlap_stats,
    parse_sa_tag,
    revcomp,
)
from tests.test_integration import (
    BAM_PATH_PAIRED_V3_0_0,
    BAM_PATH_V5_3_2,
    BED_PATH_V3_0_0,
    BED_PATH_V5_3_2,
    create_args,
)

CHROM = "MN908947.3"
HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": CHROM, "LN": 29903}]})

Primer = namedtuple(
    "Primer", ["chrom", "amplicon_number", "pool", "start", "end"], defaults=(0, 0)
)


def make_segment(start, cigar, seq=None, name="read", reverse=False, mapq=60):
    seg = pysam.AlignedSegment(HEADER)
    seg.query_name = name
    seg.reference_id = 0
    seg.reference_start = start
    seg.cigarstring = cigar
    seg.mapping_quality = mapq
    seg.is_reverse = reverse
    if seq is None:
        seg_len = seg.infer_read_length()
        seq = "".join(random.Random(start).choice("ACGT") for _ in range(seg_len))
    seg.query_sequence = seq
    return seg


def v3_qc(min_mapq=20, padding=5):
    scheme = Scheme.from_file(str(BED_PATH_V3_0_0))
    bedlines = list(scheme.bedlines)
    amplicons = create_amplicons(merge_primers(scheme.bedlines))
    lookup = create_primer_lookup([(CHROM, 29903)], amplicons, padding)
    return QC(amplicons, lookup, bedlines, min_mapq=min_mapq, padding=padding)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def records(path):
    return [r.to_string() for r in pysam.AlignmentFile(str(path), "r")]


class TestFragmentSpans(unittest.TestCase):
    def test_single_end(self):
        self.assertEqual(fragment_spans((make_segment(100, "50M"),)), [(100, 150)])

    def test_overlapping_mates_merged(self):
        segs = (make_segment(100, "50M"), make_segment(130, "50M"))
        self.assertEqual(fragment_spans(segs), [(100, 180)])

    def test_abutting_mates_merged(self):
        segs = (make_segment(150, "50M"), make_segment(100, "50M"))
        self.assertEqual(fragment_spans(segs), [(100, 200)])

    def test_disjoint_mates(self):
        segs = (make_segment(100, "50M"), make_segment(200, "50M"))
        self.assertEqual(fragment_spans(segs), [(100, 150), (200, 250)])


class TestBugFixes(unittest.TestCase):
    def test_pairing_ignores_supplementary(self):
        """A proper-pair flagged supplementary record must not be paired in place of the real mate"""
        read1 = make_segment(100, "50M", name="q")
        read1.flag = 1 + 2 + 64
        supp = make_segment(500, "50M", name="q")
        supp.flag = 1 + 2 + 64 + 2048
        read2 = make_segment(130, "50M", name="q", reverse=True)
        read2.flag = 1 + 2 + 128 + 16
        pairs = list(read_pair_generator([read1, supp, read2]))
        self.assertEqual(len(pairs), 1)
        self.assertIs(pairs[0][0], read1)
        self.assertIs(pairs[0][1], read2)

    def test_trim_invalid_cigar_raises(self):
        """trim must raise rather than return a segment whose position has moved but CIGAR has not"""
        seg = make_segment(100, "10D50M")
        with self.assertRaises(ValueError):
            trim(seg, 105, False)
        self.assertEqual(seg.cigarstring, "10D50M")

    def test_trim_removes_trailing_deletion(self):
        """A deletion left next to the right primer softmask is removed, as for the left primer"""
        seg = make_segment(100, "50M1D22M")
        trim(seg, 151, True)
        self.assertEqual(seg.cigarstring, "50M22S")
        self.assertEqual(seg.reference_end, 150)

        left = make_segment(100, "22M1D50M")
        trim(left, 123, False)
        self.assertEqual(left.cigarstring, "22S50M")
        self.assertEqual(left.reference_start, 123)

    def test_trim_extended_cigar_ops(self):
        """=/X CIGAR operations are trimmed the same as M"""
        m_seg = make_segment(100, "50M")
        x_seg = make_segment(100, "20=1X29=")
        for seg in (m_seg, x_seg):
            trim(seg, 110, False)
            trim(seg, 140, True)
        self.assertEqual(m_seg.reference_start, x_seg.reference_start)
        self.assertEqual(m_seg.reference_end, x_seg.reference_end)
        self.assertEqual(m_seg.query_alignment_start, x_seg.query_alignment_start)
        self.assertEqual(m_seg.query_alignment_end, x_seg.query_alignment_end)
        # trim re-adds the partially consumed operation as M
        self.assertEqual(x_seg.cigarstring, "10S10M1X19M10S")


class TestClassifyPrimerPair(unittest.TestCase):
    def test_categories(self):
        cases = [
            (Primer(CHROM, 5, 1), Primer(CHROM, 5, 1), "expected"),
            (Primer(CHROM, 6, 2), Primer(CHROM, 5, 1), "overlap_product"),
            (Primer(CHROM, 5, 1), Primer(CHROM, 7, 1), "long_same_pool"),
            (Primer(CHROM, 5, 1), Primer(CHROM, 6, 2), "long_cross_pool"),
            (Primer(CHROM, 5, 1), Primer(CHROM, 3, 1), "other"),
            (Primer(CHROM, 5, 1), Primer("other", 5, 1), "cross_chrom"),
        ]
        for p1, p2, expected in cases:
            self.assertEqual(classify_primer_pair(p1, p2), expected)

    def test_overlap_product_requires_primer_ends(self):
        """An overlap product must start and end within the threshold (default 5) of the primer 5' ends"""
        left = Primer(CHROM, 6, 2, start=1000, end=1022)
        right = Primer(CHROM, 5, 1, start=1080, end=1100)
        cases = [
            ((1000, 1100), "overlap_product"),
            ((1005, 1095), "overlap_product"),
            ((995, 1105), "overlap_product"),
            ((1006, 1100), "partial_in_overlap"),
            ((1000, 1094), "partial_in_overlap"),
        ]
        for (start, end), expected in cases:
            self.assertEqual(
                classify_primer_pair(left, right, start, end), expected, (start, end)
            )

    def test_long_products_not_end_checked(self):
        """Long products cannot come from a single amplicon so stay classified wherever their ends are"""
        left = Primer(CHROM, 5, 1, start=1000, end=1022)
        right = Primer(CHROM, 7, 1, start=1900, end=1922)
        self.assertEqual(
            classify_primer_pair(left, right, 1100, 1800), "long_same_pool"
        )


class TestReadPairFates(unittest.TestCase):
    def mates(self, name, flag1, flag2):
        read1 = make_segment(100, "50M", name=name)
        read1.flag = flag1
        read2 = make_segment(130, "50M", name=name, reverse=True)
        read2.flag = flag2
        return read1, read2

    def test_unpaired_records_recorded(self):
        qc = v3_qc()
        proper = self.mates("proper", 1 + 2 + 64, 1 + 2 + 128 + 16)
        not_proper = self.mates("not_proper", 1 + 64, 1 + 128 + 16)
        # read2 unmapped, read1 flags its mate as unmapped
        unmapped = self.mates("unmapped", 1 + 8 + 64, 1 + 4 + 128)
        orphan = self.mates("orphan", 1 + 2 + 64, 1 + 2 + 128 + 16)[0]
        orphan_mate_unmapped = self.mates("orphan_unmapped", 1 + 8 + 64, 0)[0]

        records = [*proper, *not_proper, *unmapped, orphan, orphan_mate_unmapped]
        pairs = list(read_pair_generator(records, qc=qc))

        self.assertEqual(pairs, [proper])
        self.assertEqual(
            dict(qc.unassigned),
            {"not_proper_pair": 1, "unmapped": 2, "missing_mate": 1},
        )

    def test_supplementary_and_secondary_not_counted(self):
        qc = v3_qc()
        supp = make_segment(500, "50M", name="proper")
        supp.flag = 1 + 2 + 64 + 2048
        secondary = make_segment(700, "50M", name="proper")
        secondary.flag = 1 + 2 + 128 + 256
        proper = self.mates("proper", 1 + 2 + 64, 1 + 2 + 128 + 16)
        pairs = list(
            read_pair_generator([proper[0], supp, secondary, proper[1]], qc=qc)
        )
        self.assertEqual(pairs, [proper])
        self.assertEqual(sum(qc.unassigned.values()), 0)


class TestSplitAlignments(unittest.TestCase):
    def test_parse_sa_tag(self):
        parts = parse_sa_tag("chr1,101,+,10S50M40S,60,2;chr2,201,-,10H5M2D45M40H,7,0;")
        self.assertEqual(len(parts), 2)
        self.assertEqual(parts[0], ("chr1", 100, 150, "+", 60, 10, 60))
        # minus strand query coordinates are converted to read orientation
        self.assertEqual(parts[1], ("chr2", 200, 252, "-", 7, 40, 90))

    def test_junctions_positions_and_orientation(self):
        primary = parse_sa_tag("c,1,+,100M50S,60,0")[0]
        collinear = parse_sa_tag("c,501,+,100S50M,60,0")[0]
        self.assertEqual(junctions([collinear, primary]), [(primary, collinear)])
        # last aligned base before the break (100) and first after it (501), 1-based
        self.assertEqual(breakpoint_positions(primary, collinear), (100, 501))
        self.assertEqual(join_orientation(primary, collinear), "collinear")
        backwards = parse_sa_tag("c,1001,+,100S50M,60,0")[0]._replace(ref_start=-50)
        self.assertEqual(join_orientation(primary, backwards), "non_collinear")
        inverted = parse_sa_tag("c,501,-,50M100S,60,0")[0]
        self.assertEqual(join_orientation(primary, inverted), "inverted")
        # on the minus strand the part's read-order end is its lowest reference base
        self.assertEqual(breakpoint_positions(primary, inverted), (100, 550))

    def classify(self, sa):
        qc = v3_qc()
        seg = make_segment(30, "380M300S")
        seg.set_tag("SA", sa)
        qc.split((seg,))
        self.assertEqual(sum(qc.split_fragments.values()), 1)
        return next(iter(qc.split_fragments))

    def test_amplicon_chimera(self):
        # amplicon 1 (30-410) joined to amplicon 3 (starts 642)
        self.assertEqual(
            self.classify(f"{CHROM},643,+,380S300M,60,0;"), "amplicon_chimera"
        )

    def test_within_amplicon(self):
        # 200 is inside the amplicon 1 insert, away from every amplicon end
        self.assertEqual(
            self.classify(f"{CHROM},201,+,380S300M,60,0;"), "within_amplicon"
        )

    def test_low_mapq(self):
        self.assertEqual(self.classify(f"{CHROM},643,+,380S300M,5,0;"), "low_mapq")

    def test_no_sa_tag(self):
        qc = v3_qc()
        qc.split((make_segment(30, "380M"),))
        self.assertEqual(qc.split_checked, 1)
        self.assertEqual(sum(qc.split_fragments.values()), 0)

    def test_join_in_both_mates_counted_once(self):
        """Overlapping mates of one pair can both cross the same join, which is one event"""
        qc = v3_qc()
        mate1 = make_segment(30, "380M300S")
        mate1.set_tag("SA", f"{CHROM},643,+,380S300M,60,0;")
        # the reverse mate reads the same molecule from the other end: amplicon 3 part first in read order
        mate2 = make_segment(642, "300S300M", reverse=True)
        mate2.set_tag("SA", f"{CHROM},31,-,380M300S,60,0;")
        qc.split((mate1, mate2))
        self.assertEqual(qc.split_checked, 1)
        self.assertEqual(dict(qc.split_fragments), {"amplicon_chimera": 1})
        self.assertEqual(list(qc.joins.values()), [1])


class TestSameDirection(unittest.TestCase):
    def test_left_left_detected(self):
        qc = v3_qc()
        left_3 = next(
            bl
            for bl in qc.left_end_kmers.values()
            if bl.primername == "nCoV-2019_3_LEFT_1"
        )
        # read starts at LEFT_1, ends well inside amplicon 1 carrying LEFT_3 reverse complemented
        tail = revcomp(left_3.sequence[-12:])
        body = "".join(random.Random(1).choice("ACGT") for _ in range(300 - len(tail)))
        seg = make_segment(30, "300M", seq=body + tail)
        p1 = qc._find_primer(qc.lookup, 30, "+", CHROM)
        p2 = qc._find_primer(qc.lookup, seg.reference_end, "-", CHROM)
        qc.same_direction(seg, seg, p1, p2)
        self.assertEqual(
            [(a.primername, b.primername, c) for a, b, c in qc.pairs],
            [("nCoV-2019_1_LEFT_1", "nCoV-2019_3_LEFT_1", "two_left_primers")],
        )

    def test_explained_ends_not_scanned(self):
        qc = v3_qc()
        seg = make_segment(30, "380M")
        p1 = qc._find_primer(qc.lookup, 30, "+", CHROM)
        p2 = qc._find_primer(qc.lookup, seg.reference_end, "-", CHROM)
        qc.same_direction(seg, seg, p1, p2)
        self.assertEqual(len(qc.pairs), 0)


class TestOverlapCounting(unittest.TestCase):
    def test_count_bases_in_window(self):
        seq = "TTTTT" + "ACGTACGTAC" + "GGGGG" + "T" + "CCCCCCCCCC"
        seg = make_segment(100, "5S10M2D5M1I10M", seq=seq)
        counts = np.zeros((30, 7), dtype=int)
        covered = count_bases_in_window(seg, 100, 130, counts)
        # soft clip excluded, so the first window base is the first M base
        self.assertEqual(list(counts[0]), [1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(list(counts[10]), [0, 0, 0, 0, 0, 1, 0])  # deletion
        self.assertEqual(counts[16, 6], 1)  # insertion after the 5M block
        self.assertEqual(int(covered.sum()), 27)
        self.assertEqual(int(counts[:, :6].sum()), 27)

        # the other mate of the same template does not count the same positions again
        count_bases_in_window(seg, 100, 130, counts, seen=covered)
        self.assertEqual(int(counts[:, :6].sum()), 27)

    def test_overlap_stats(self):
        left = np.zeros((2, 7), dtype=int)
        right = np.zeros((2, 7), dtype=int)
        left[0, 0], right[0, 0] = 20, 20  # concordant A
        left[1, 0], left[1, 3], right[1, 3] = (
            10,
            10,
            5,
        )  # 50% A/T vs 100% T, right depth too low
        _, _, af_diff, evaluable, left_major, right_major = overlap_stats(left, right)
        self.assertEqual(list(evaluable), [True, False])
        self.assertAlmostEqual(af_diff[1], 0.5)
        self.assertEqual(af_diff[0], 0)
        self.assertEqual(list(left_major), ["A", "A"])
        self.assertEqual(right_major[1], "T")


class TestBaseFrequencies(unittest.TestCase):
    def test_most_common_first(self):
        counts = np.array([5, 15, 0, 0, 0, 0, 3])  # A C G T N DEL INS
        self.assertEqual(base_frequencies(counts), "C=0.75,A=0.25")

    def test_deletions_counted_and_n_ignored(self):
        counts = np.array([0, 0, 0, 9, 4, 1, 0])
        self.assertEqual(base_frequencies(counts), "T=0.90,DEL=0.10")

    def test_rare_bases_left_out(self):
        counts = np.array([999, 0, 1, 0, 0, 0, 0])
        self.assertEqual(base_frequencies(counts), "A=1.00")

    def test_details_column(self):
        qc = v3_qc()
        window = qc.windows[0]
        window["left"][0, 0], window["left"][0, 3] = 18, 2  # A 0.9, T 0.1
        window["right"][0, 3], window["right"][0, 0] = 15, 5  # T 0.75, A 0.25
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-freq") as tempdir:
            path = f"{tempdir}/overlaps.tsv"
            with contextlib.redirect_stderr(io.StringIO()):
                console.reset()
                qc._write_overlaps(path)
            row = read_tsv(path)[0]
        self.assertEqual(row["majority_base_differences"], "1")
        self.assertEqual(
            row["majority_base_difference_details"],
            f"{int(row['start'])}:A=0.90,T=0.10/T=0.75,A=0.25",
        )


class TestContaminationWarning(unittest.TestCase):
    def write_overlaps(self, qc):
        stderr = io.StringIO()
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-warn") as tempdir:
            with contextlib.redirect_stderr(stderr):
                console.reset()
                qc._write_overlaps(f"{tempdir}/overlaps.tsv")
                # the warning is held back and printed under the Warnings heading
                console.flush_warnings()
        return stderr.getvalue()

    def test_warning_on_major_allele_mismatch(self):
        qc = v3_qc()
        window = qc.windows[0]
        window["left"][0, 0] = 20  # A
        window["right"][0, 1] = 20  # C
        out = self.write_overlaps(qc)
        self.assertIn("WARNING: majority base differs", out)
        self.assertIn("which can indicate contamination", out)
        self.assertIn("1/2", out)

    def test_no_warning_when_concordant(self):
        qc = v3_qc()
        window = qc.windows[0]
        window["left"][0, 0] = 20
        window["right"][0, 0] = 20
        self.assertNotIn("WARNING", self.write_overlaps(qc))


class TestQCIntegration(unittest.TestCase):
    def run_go(self, tempdir, name, **kwargs):
        out = pathlib.Path(tempdir) / f"{name}.sam"
        args = create_args(output=out.absolute(), **kwargs)
        go(args)
        return out

    def test_output_unchanged_by_qc(self):
        for kwargs in (
            {"normalise": 0},
            {"normalise": 20},
            {
                "normalise": 0,
                "bedfile": BED_PATH_V3_0_0,
                "samfile": BAM_PATH_PAIRED_V3_0_0,
            },
            {
                "normalise": 200,
                "bedfile": BED_PATH_V3_0_0,
                "samfile": BAM_PATH_PAIRED_V3_0_0,
            },
        ):
            with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-same") as tempdir:
                plain = self.run_go(tempdir, "plain", no_read_groups=True, **kwargs)
                with_qc = self.run_go(
                    tempdir,
                    "qc",
                    no_read_groups=True,
                    qc_prefix=f"{tempdir}/qc",
                    **kwargs,
                )
                self.assertEqual(records(plain), records(with_qc), kwargs)

    def assert_fragments_add_up(self, rows, reasons):
        for r in rows:
            removed = sum(int(r[f"removed_{reason}"]) for reason in reasons)
            self.assertEqual(int(r["kept"]) + removed, int(r["fragments"]), r)

    def test_fragment_accounting_se(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-se") as tempdir:
            for normalise in (0, 20):
                prefix = f"{tempdir}/n{normalise}"
                report = pathlib.Path(tempdir) / f"n{normalise}.report.tsv"
                out = self.run_go(
                    tempdir,
                    f"n{normalise}",
                    normalise=normalise,
                    qc_prefix=prefix,
                    report=report,
                )
                (summary,) = read_tsv(f"{prefix}.summary.tsv")
                amplicons = read_tsv(f"{prefix}.amplicons.tsv")
                self.assert_fragments_add_up([summary], REASONS)
                self.assert_fragments_add_up(amplicons, AMPLICON_REASONS)

                # each read counted once, extra alignments of a read are not
                self.assertEqual(
                    int(summary["fragments"]),
                    sum(
                        1
                        for r in pysam.AlignmentFile(str(BAM_PATH_V5_3_2))
                        if not (r.is_supplementary or r.is_secondary)
                    ),
                )
                self.assertEqual(int(summary["kept"]), len(records(out)))
                self.assertEqual(
                    sum(int(r["kept"]) for r in amplicons), int(summary["kept"])
                )
                self.assertEqual(
                    int(summary["removed_by_normalise"]) > 0, normalise > 0
                )

                scheme_amplicons = create_amplicons(
                    merge_primers(Scheme.from_file(create_args().bedfile).bedlines)
                )
                self.assertEqual(
                    [r["amplicon"] for r in amplicons],
                    [str(a.amplicon_number) for a in scheme_amplicons],
                )

                mismatched = sum(
                    1 for r in read_tsv(report) if r["CorrectlyPaired"] == "False"
                )
                self.assertEqual(
                    int(summary["removed_mismatched_primers"])
                    + int(summary["removed_ambiguous_amplicon"]),
                    mismatched,
                )

    def test_fragment_accounting_pe(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-pe") as tempdir:
            for normalise in (0, 200):
                prefix = f"{tempdir}/n{normalise}"
                out = self.run_go(
                    tempdir,
                    f"n{normalise}",
                    normalise=normalise,
                    qc_prefix=prefix,
                    bedfile=BED_PATH_V3_0_0,
                    samfile=BAM_PATH_PAIRED_V3_0_0,
                )
                (summary,) = read_tsv(f"{prefix}.summary.tsv")
                self.assert_fragments_add_up([summary], REASONS)
                self.assert_fragments_add_up(
                    read_tsv(f"{prefix}.amplicons.tsv"), AMPLICON_REASONS
                )
                self.assertEqual(int(summary["kept"]) * 2, len(records(out)))
                # every fragment (read name) is counted exactly once, paired or not
                fragments = {
                    r.query_name
                    for r in pysam.AlignmentFile(str(BAM_PATH_PAIRED_V3_0_0))
                    if not (r.is_supplementary or r.is_secondary)
                }
                self.assertEqual(int(summary["fragments"]), len(fragments))

    def test_primer_pairs(self):
        scheme = Scheme.from_file(create_args().bedfile)
        # CorrectlyPaired ignores chrom, so it only equals "expected" for single chrom schemes
        self.assertEqual(len({bl.chrom for bl in scheme.bedlines}), 1)

        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-pairs") as tempdir:
            counts = {}
            for allow in (False, True):
                prefix = f"{tempdir}/allow{allow}"
                report = pathlib.Path(tempdir) / f"allow{allow}.report.tsv"
                self.run_go(
                    tempdir,
                    f"allow{allow}",
                    allow_incorrect_pairs=allow,
                    qc_prefix=prefix,
                    report=report,
                )
                rows = read_tsv(f"{prefix}.primer_pairs.tsv")
                counts[allow] = sorted(
                    (r["left_primer"], r["right_primer"], r["type"], r["fragments"])
                    for r in rows
                )
                coordinate_rows = [
                    r
                    for r in rows
                    if r["type"] not in ("two_left_primers", "two_right_primers")
                ]
                report_rows = read_tsv(report)
                self.assertEqual(
                    sum(int(r["fragments"]) for r in coordinate_rows), len(report_rows)
                )
                self.assertEqual(
                    sum(int(r["fragments"]) for r in rows if r["type"] == "expected"),
                    sum(1 for r in report_rows if r["CorrectlyPaired"] == "True"),
                )
            self.assertEqual(counts[False], counts[True])

    def test_split_reads_se(self):
        expected = sum(
            1
            for r in pysam.AlignmentFile(str(BAM_PATH_V5_3_2))
            if not (r.is_unmapped or r.is_supplementary or r.is_secondary)
            and r.mapping_quality >= 20
            and r.has_tag("SA")
        )
        self.assertGreater(expected, 0)
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-split") as tempdir:
            prefix = f"{tempdir}/qc"
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                self.run_go(tempdir, "qc", qc_prefix=prefix)
            self.assertIn(f"QC: split reads: {expected} of ", stderr.getvalue())
            rows = read_tsv(f"{prefix}.split_reads.tsv")
            # every split read has at least one join
            self.assertGreaterEqual(sum(int(r["fragments"]) for r in rows), expected)
            self.assertTrue(
                {r["type"] for r in rows}
                <= {"amplicon_chimera", "within_amplicon", "low_mapq"}
            )
            for r in rows:
                # lower coordinate side first, positions are 1-based window starts
                self.assertLessEqual(
                    (r["chrom"], int(r["pos"])),
                    (r["joined_chrom"], int(r["joined_pos"])),
                )
                self.assertEqual(int(r["pos"]) % 10, 1)

    def test_overlaps(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-overlap") as tempdir:
            prefix = f"{tempdir}/qc"
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                self.run_go(
                    tempdir,
                    "qc",
                    qc_prefix=prefix,
                    bedfile=BED_PATH_V3_0_0,
                    samfile=BAM_PATH_PAIRED_V3_0_0,
                )
            overlaps = read_tsv(f"{prefix}.overlaps.tsv")
            compared = [r for r in overlaps if int(r["compared_positions"])]
            self.assertTrue(compared)
            # the warning is printed exactly when an overlap has a majority base difference
            self.assertEqual(
                "WARNING: majority base differs" in stderr.getvalue(),
                any(int(r["majority_base_differences"]) for r in compared),
            )
            for r in compared:
                details = r["majority_base_difference_details"]
                self.assertEqual(
                    len(details.split(";")) if details else 0,
                    int(r["majority_base_differences"]),
                )
            for r in overlaps:
                if not int(r["compared_positions"]):
                    self.assertEqual(r["majority_base_differences"], "NA")
            first = overlaps[0]
            self.assertEqual((first["amplicon"], first["next_amplicon"]), ("1", "2"))
            # [342, 385) 0-based is 343-385 1-based inclusive
            self.assertEqual((first["start"], first["end"]), ("343", "385"))


class TestFragmentNormalisation(unittest.TestCase):
    def amp_depths(self, tempdir, normalise):
        path = pathlib.Path(tempdir) / f"amp_{normalise}.tsv"
        args = create_args(
            output=(pathlib.Path(tempdir) / f"out_{normalise}.bam").absolute(),
            bedfile=BED_PATH_V3_0_0,
            samfile=BAM_PATH_PAIRED_V3_0_0,
            amp_depth_report=path,
            normalise=normalise,
        )
        go(args)
        return {r["amplicon"]: float(r["mean_fragment_depth"]) for r in read_tsv(path)}

    def test_paired_normalisation_counts_fragments(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-fragments") as tempdir:
            depths = [d for d in self.amp_depths(tempdir, 100).values() if d > 0]
            self.assertTrue(70 <= statistics.median(depths) <= 130)

    def test_paired_amp_depth_counts_both_mates(self):
        """Without normalisation the paired amp depth report counts whole fragments (it used to count one mate)"""
        with tempfile.TemporaryDirectory(dir="tests", suffix="-amp-depth") as tempdir:
            no_norm = self.amp_depths(tempdir, 0)
            keep_all = self.amp_depths(tempdir, 1_000_000)
            for amplicon, depth in no_norm.items():
                self.assertAlmostEqual(depth, keep_all[amplicon], places=6)


def primary_records(path):
    return [
        r
        for r in pysam.AlignmentFile(str(path))
        if not (r.is_supplementary or r.is_secondary)
    ]


class FakeWriter:
    def __init__(self):
        self.records = []

    def write(self, record):
        self.records.append(record)


class TestRejectedOutput(unittest.TestCase):
    RUNS = (
        {"normalise": 0},
        {"normalise": 20},
        {"normalise": 0, "bedfile": BED_PATH_V3_0_0, "samfile": BAM_PATH_PAIRED_V3_0_0},
        {
            "normalise": 200,
            "bedfile": BED_PATH_V3_0_0,
            "samfile": BAM_PATH_PAIRED_V3_0_0,
        },
    )

    def run_go(self, tempdir, name, **kwargs):
        out = pathlib.Path(tempdir) / f"{name}.sam"
        with contextlib.redirect_stderr(io.StringIO()):
            go(create_args(output=out.absolute(), **kwargs))
        return out

    def test_every_fragment_kept_or_rejected(self):
        for kwargs in self.RUNS:
            samfile = kwargs.get("samfile", create_args().samfile)
            with tempfile.TemporaryDirectory(
                dir="tests", suffix="-rejected"
            ) as tempdir:
                rejected_path = f"{tempdir}/rejected.bam"
                out = self.run_go(
                    tempdir,
                    "out",
                    qc_prefix=f"{tempdir}/qc",
                    rejected_output=rejected_path,
                    **kwargs,
                )
                inputs = primary_records(samfile)
                kept = list(pysam.AlignmentFile(str(out)))
                rejected = list(pysam.AlignmentFile(rejected_path))

                # every primary input record is in exactly one of the two files
                key = lambda r: (r.query_name, r.flag)  # noqa: E731
                self.assertEqual(
                    Counter(map(key, inputs)),
                    Counter(map(key, kept)) + Counter(map(key, rejected)),
                    kwargs,
                )
                self.assertFalse(
                    {r.query_name for r in kept} & {r.query_name for r in rejected}
                )
                self.assertFalse(
                    any(r.is_supplementary or r.is_secondary for r in rejected)
                )

                # one fragment per read name, counted under the same reasons as the QC summary
                reasons = Counter(
                    dict((r.query_name, r.get_tag("ZR")) for r in rejected).values()
                )
                (summary,) = read_tsv(f"{tempdir}/qc.summary.tsv")
                self.assertEqual(
                    {r: reasons[r] for r in REASONS},
                    {r: int(summary[f"removed_{r}"]) for r in REASONS},
                    kwargs,
                )

                # written as read in, apart from the reason tag
                originals = {key(r): r.to_string() for r in inputs}
                for r in rejected:
                    r.set_tag("ZR", None)
                    self.assertEqual(r.to_string(), originals[key(r)])

                # the same reasons without --qc-prefix
                self.run_go(
                    tempdir,
                    "no_qc",
                    rejected_output=f"{tempdir}/no_qc.rejected.sam",
                    **kwargs,
                )
                self.assertEqual(
                    Counter(
                        r.get_tag("ZR")
                        for r in pysam.AlignmentFile(f"{tempdir}/no_qc.rejected.sam")
                    ),
                    Counter(
                        r.get_tag("ZR") for r in pysam.AlignmentFile(rejected_path)
                    ),
                )

    def test_output_unchanged(self):
        for kwargs in self.RUNS:
            with tempfile.TemporaryDirectory(
                dir="tests", suffix="-rejected"
            ) as tempdir:
                plain = self.run_go(tempdir, "plain", **kwargs)
                with_rejected = self.run_go(
                    tempdir,
                    "with",
                    rejected_output=f"{tempdir}/rejected.bam",
                    **kwargs,
                )
                self.assertEqual(records(plain), records(with_rejected), kwargs)

    def test_header(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-rejected") as tempdir:
            path = f"{tempdir}/rejected.sam"
            self.run_go(
                tempdir,
                "out",
                rejected_output=path,
                bedfile=BED_PATH_V3_0_0,
                samfile=BAM_PATH_PAIRED_V3_0_0,
            )
            header = pysam.AlignmentFile(path).header.to_dict()
            self.assertEqual(header["HD"]["SO"], "unsorted")
            self.assertEqual(header["PG"][-1]["ID"], "align_trim")

    def test_read_pair_generator(self):
        def mate(name, start, read1, proper=True, unmapped=False, supplementary=False):
            seg = make_segment(start, "50M", name=name, reverse=not read1)
            seg.is_paired = True
            seg.is_read1, seg.is_read2 = read1, not read1
            seg.is_proper_pair = proper
            seg.is_unmapped = unmapped
            seg.is_supplementary = supplementary
            return seg

        reads = [
            mate("proper", 100, True),
            mate("proper", 200, False),
            mate("proper", 5000, True, supplementary=True),
            mate("improper", 300, True, proper=False),
            mate("improper", 9000, False, proper=False),
            mate("unmapped", 400, True, proper=False),
            mate("unmapped", 400, False, proper=False, unmapped=True),
            mate("orphan", 500, True),
        ]
        writer = FakeWriter()
        pairs = list(read_pair_generator(reads, rejected=writer))
        self.assertEqual([p[0].query_name for p in pairs], ["proper"])
        self.assertEqual(
            [(r.query_name, r.get_tag("ZR")) for r in writer.records],
            [
                ("improper", "not_proper_pair"),
                ("improper", "not_proper_pair"),
                ("unmapped", "unmapped"),
                ("unmapped", "unmapped"),
                ("orphan", "missing_mate"),
            ],
        )

    def test_invalid_path(self):
        for path in ("-", "rejected.txt"):
            with tempfile.TemporaryDirectory(
                dir="tests", suffix="-rejected"
            ) as tempdir:
                with self.assertRaises(SystemExit) as cm:
                    self.run_go(tempdir, "out", rejected_output=pathlib.Path(path))
                self.assertEqual(cm.exception.code, 1)


def v532(padding=5, ref_len=29903):
    scheme = Scheme.from_file(str(BED_PATH_V5_3_2))
    bedlines = list(scheme.bedlines)
    amplicons = create_amplicons(merge_primers(scheme.bedlines))
    lookup = create_primer_lookup([(CHROM, ref_len)], amplicons, padding)
    return amplicons, lookup, QC(amplicons, lookup, bedlines, 20, padding)


class TestRightEndLookup(unittest.TestCase):
    """A right end within the threshold of its primer's 5' end matches it, as on the left"""

    def test_threshold_is_inclusive_on_both_sides(self):
        amplicons, lookup, _ = v532()
        amp1 = amplicons[0]  # LEFT 47-78, RIGHT 419-447
        right = find_primer_with_lookup(lookup, amp1.amplicon_end + 5, "-", CHROM)
        left = find_primer_with_lookup(lookup, amp1.amplicon_start - 5, "+", CHROM)
        self.assertEqual((left.amplicon_number, right.amplicon_number), (1, 1))
        beyond = find_primer_with_lookup(lookup, amp1.amplicon_end + 6, "-", CHROM)
        self.assertNotEqual(beyond.amplicon_number, 1)

    def test_read_ending_at_the_reference_end(self):
        amplicons, _, _ = v532()
        last = amplicons[-1]
        _, lookup, _ = v532(ref_len=last.amplicon_end)
        p = find_primer_with_lookup(lookup, last.amplicon_end, "-", CHROM)
        self.assertEqual(p.amplicon_number, last.amplicon_number)


class TestAmbiguousAmplicon(unittest.TestCase):
    """Fragments inside two neighbouring amplicons are not assigned to either"""

    def handle(self, seg, **kwargs):
        amplicons, lookup, qc = v532()
        rejected = FakeWriter()
        amp_depths = {
            CHROM: {
                a.amplicon_number: np.zeros(
                    a.amplicon_end - a.amplicon_start, dtype=int
                )
                for a in amplicons
            }
        }
        result = handle_segments(
            segment=seg,
            lookup=lookup,
            args=create_args(**kwargs),
            min_mapq=20,
            outfile_writer=FakeWriter(),
            amp_depths=amp_depths,
            qc=qc,
            rejected=rejected,
            originals=(copy.copy(seg),),
        )
        return result, qc, rejected

    def assert_ambiguous(self, seg, **kwargs):
        result, qc, rejected = self.handle(seg, **kwargs)
        self.assertFalse(result)
        self.assertEqual(qc.unassigned["ambiguous_amplicon"], 1)
        self.assertEqual(sum(sum(counts.values()) for counts in qc.fates.values()), 0)
        self.assertEqual(qc.ambiguous, Counter({(CHROM, 1): 1}))
        self.assertEqual(
            [r.get_tag("ZR") for r in rejected.records], ["ambiguous_amplicon"]
        )

    def test_partial_read_in_overlap(self):
        # 354-447 (1-based): ends at RIGHT_1, starts inside LEFT_2, so lies inside amplicons 1 and 2
        for allow in (False, True):
            seg = make_segment(353, "94M26S", reverse=True)
            self.assert_ambiguous(seg, allow_incorrect_pairs=allow)

    def test_overlap_product(self):
        # both ends at LEFT_2 / RIGHT_1
        for allow in (False, True):
            self.assert_ambiguous(
                make_segment(344, "103M"), allow_incorrect_pairs=allow
            )

    def test_counted_per_overlap(self):
        _, _, qc = v532()
        qc.ambiguous[(CHROM, 1)] = 3
        # 67/68 only overlap within their primers, so have no comparison region
        qc.ambiguous[(CHROM, 67)] = 2
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-amb") as tempdir:
            path = f"{tempdir}/overlaps.tsv"
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                qc._write_overlaps(path)
            rows = {(r["amplicon"], r["next_amplicon"]): r for r in read_tsv(path)}
        self.assertEqual(rows[("1", "2")]["ambiguous_fragments"], "3")
        self.assertEqual(rows[("2", "3")]["ambiguous_fragments"], "0")
        self.assertEqual(rows[("67", "68")]["ambiguous_fragments"], "2")
        self.assertEqual(rows[("67", "68")]["start"], "NA")
        self.assertIn("5 fragments lie inside two amplicons", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()

import csv
import pathlib
import random
import statistics
import tempfile
import unittest
from collections import namedtuple

import numpy as np
import pysam
from primalbedtools.amplicons import create_amplicons
from primalbedtools.bedfiles import merge_primers
from primalbedtools.scheme import Scheme

from align_trim.main import (
    create_primer_lookup,
    fragment_spans,
    go,
    read_pair_generator,
    trim,
)
from align_trim.qc import (
    QC,
    classify_primer_pair,
    count_bases_in_window,
    junctions,
    overlap_stats,
    parse_sa_tag,
    revcomp,
    split_orientation,
)
from tests.test_integration import (
    BAM_PATH_PAIRED_V3_0_0,
    BAM_PATH_V5_3_2,
    BED_PATH_V3_0_0,
    create_args,
)

CHROM = "MN908947.3"
HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": CHROM, "LN": 29903}]})

Primer = namedtuple("Primer", ["chrom", "amplicon_number", "pool"])


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


def v3_qc(min_mapq=20, padding=35):
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
            (Primer(CHROM, 5, 1), Primer(CHROM, 5, 1), "correct"),
            (Primer(CHROM, 6, 2), Primer(CHROM, 5, 1), "short_overlap"),
            (Primer(CHROM, 5, 1), Primer(CHROM, 7, 1), "long_same_pool"),
            (Primer(CHROM, 5, 1), Primer(CHROM, 6, 2), "long_cross_pool"),
            (Primer(CHROM, 5, 1), Primer(CHROM, 3, 1), "other"),
            (Primer(CHROM, 5, 1), Primer("other", 5, 1), "cross_chrom"),
        ]
        for p1, p2, expected in cases:
            self.assertEqual(classify_primer_pair(p1, p2), expected)


class TestSplitAlignments(unittest.TestCase):
    def test_parse_sa_tag(self):
        parts = parse_sa_tag("chr1,101,+,10S50M40S,60,2;chr2,201,-,10H5M2D45M40H,7,0;")
        self.assertEqual(len(parts), 2)
        self.assertEqual(parts[0], ("chr1", 100, 150, "+", 60, 10, 60))
        # minus strand query coordinates are converted to read orientation
        self.assertEqual(parts[1], ("chr2", 200, 252, "-", 7, 40, 90))

    def test_junctions_and_orientation(self):
        primary = parse_sa_tag("c,1,+,100M50S,60,0")[0]
        collinear = parse_sa_tag("c,501,+,100S50M,60,0")[0]
        self.assertEqual(junctions([collinear, primary]), [(("c", 100), ("c", 500))])
        self.assertEqual(split_orientation(primary, [primary, collinear]), "collinear")
        backwards = parse_sa_tag("c,1001,+,100S50M,60,0")[0]._replace(ref_start=-50)
        self.assertEqual(
            split_orientation(primary, [primary, backwards]), "non_collinear"
        )
        inverted = parse_sa_tag("c,501,-,50M100S,60,0")[0]
        self.assertEqual(split_orientation(primary, [primary, inverted]), "inverted")

    def classify(self, sa):
        qc = v3_qc()
        seg = make_segment(30, "380M300S")
        seg.set_tag("SA", sa)
        qc.split(seg)
        self.assertEqual(sum(qc.split_classes.values()), 1)
        return next(iter(qc.split_classes))

    def test_amplicon_junction(self):
        # amplicon 1 (30-410) joined to amplicon 3 (starts 642)
        self.assertEqual(
            self.classify(f"{CHROM},643,+,380S300M,60,0;"), "amplicon_junction"
        )

    def test_internal_breakpoint(self):
        # 200 is inside the amplicon 1 insert, away from every amplicon end
        self.assertEqual(
            self.classify(f"{CHROM},201,+,380S300M,60,0;"), "internal_breakpoint"
        )

    def test_low_confidence(self):
        self.assertEqual(
            self.classify(f"{CHROM},643,+,380S300M,5,0;"), "low_confidence"
        )

    def test_no_sa_tag(self):
        qc = v3_qc()
        qc.split(make_segment(30, "380M"))
        self.assertEqual(qc.split_reads, 1)
        self.assertEqual(sum(qc.split_classes.values()), 0)


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
            [("nCoV-2019_1_LEFT_1", "nCoV-2019_3_LEFT_1", "left_left")],
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

    def test_read_fate_se(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-fate") as tempdir:
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
                rows = read_tsv(f"{prefix}.read_fate.tsv")
                total = sum(int(r["total"]) for r in rows)
                self.assertEqual(
                    total, sum(1 for _ in pysam.AlignmentFile(str(BAM_PATH_V5_3_2)))
                )
                self.assertEqual(sum(int(r["kept"]) for r in rows), len(records(out)))
                normalised_out = sum(int(r["normalised_out"]) for r in rows)
                self.assertEqual(normalised_out > 0, normalise > 0)

                amplicons = create_amplicons(
                    merge_primers(Scheme.from_file(create_args().bedfile).bedlines)
                )
                assigned = {
                    r["amplicon"] for r in rows if r["amplicon"] != "unassigned"
                }
                self.assertEqual(assigned, {str(a.amplicon_number) for a in amplicons})

                incorrect = sum(
                    1 for r in read_tsv(report) if r["CorrectlyPaired"] == "False"
                )
                self.assertEqual(
                    sum(int(r["incorrectly_paired"]) for r in rows), incorrect
                )

    def test_read_fate_pe(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-fate-pe") as tempdir:
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
                rows = read_tsv(f"{prefix}.read_fate.tsv")
                self.assertEqual(
                    sum(int(r["kept"]) for r in rows) * 2, len(records(out))
                )

    def test_primer_pairs(self):
        scheme = Scheme.from_file(create_args().bedfile)
        # CorrectlyPaired ignores chrom, so it only equals "correct" for single chrom schemes
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
                    (r["left_primer"], r["right_primer"], r["category"], r["count"])
                    for r in rows
                )
                coordinate_rows = [
                    r for r in rows if r["category"] not in ("left_left", "right_right")
                ]
                report_rows = read_tsv(report)
                self.assertEqual(
                    sum(int(r["count"]) for r in coordinate_rows), len(report_rows)
                )
                self.assertEqual(
                    sum(int(r["count"]) for r in rows if r["category"] == "correct"),
                    sum(1 for r in report_rows if r["CorrectlyPaired"] == "True"),
                )
            self.assertEqual(counts[False], counts[True])

    def test_split_alignments_se(self):
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
            self.run_go(tempdir, "qc", qc_prefix=prefix)
            rows = read_tsv(f"{prefix}.split_alignments.tsv")
            self.assertEqual(sum(int(r["reads"]) for r in rows), expected)
            self.assertTrue(
                {r["class"] for r in rows}
                <= {"amplicon_junction", "internal_breakpoint", "low_confidence"}
            )
            self.assertTrue(read_tsv(f"{prefix}.split_breakpoints.tsv"))

    def test_overlap_windows(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-qc-overlap") as tempdir:
            prefix = f"{tempdir}/qc"
            self.run_go(
                tempdir,
                "qc",
                qc_prefix=prefix,
                bedfile=BED_PATH_V3_0_0,
                samfile=BAM_PATH_PAIRED_V3_0_0,
            )
            summary = read_tsv(f"{prefix}.overlap_summary.tsv")
            first = summary[0]
            self.assertEqual(
                (first["amplicon_left"], first["amplicon_right"]), ("1", "2")
            )
            # [342, 385) 0-based is 343-385 1-based inclusive
            self.assertEqual((first["start"], first["end"]), ("343", "385"))
            positions = read_tsv(f"{prefix}.overlap_positions.tsv")
            self.assertEqual(len(positions), sum(int(r["length"]) for r in summary))
            self.assertTrue(any(int(r["evaluable_positions"]) for r in summary))


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


if __name__ == "__main__":
    unittest.main()

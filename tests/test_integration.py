import pathlib
import sys
import unittest
import argparse
from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import merge_primers
from primalbedtools.amplicons import create_amplicons
import tempfile
import pysam
from align_trim.main import go, create_primer_lookup, find_primer_with_lookup
from collections import defaultdict

BED_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/v5.3.2.primer.bed"
BAM_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/sars-cov-2_v5.3.2.bam"
BED_PATH_V3_0_0 = pathlib.Path(__file__).parent / "test_data/v3.0.0.primer.bed"
BAM_PATH_PAIRED_V3_0_0 = (
    pathlib.Path(__file__).parent / "test_data/sars-cov-2_v3.0.0_paired.bam"
)


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def create_args(**kwargs):
    """Create a fake args object with default values matching main.py argument parser"""
    defaults = {
        "bedfile": str(BED_PATH_V5_3_2.absolute()),
        "samfile": str(BAM_PATH_V5_3_2.absolute()),
        "normalise": 0,
        "min_mapq": 20,
        "primer_match_threshold": 35,
        "report": None,
        "amp_depth_report": None,
        "no_trim_primers": False,
        "paired": False,
        "no_read_groups": False,
        "verbose": False,
        "allow_incorrect_pairs": False,
        "require_full_length": False,
        "output": None,
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


class TestIntegration(unittest.TestCase):
    def test_align_trim_trim_se_no_norm(self):
        """Tests primers are trimmed correctly"""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-trim_primers_se_no_norm"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"

            # Create args with test-specific values
            args = create_args(output=output_sam.absolute())

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)
            amplicons = create_amplicons(scheme.bedlines)
            pools = set([bl.pool for bl in scheme.bedlines])

            ref_lengths = [("MN908947.3", 29903)]
            primer_lookup = create_primer_lookup(ref_lengths, amplicons, 35)

            # Check the out sam is as expected
            for record in pysam.AlignmentFile(str(output_sam), "r"):
                # Find the left primer
                lp = find_primer_with_lookup(
                    primer_lookup, record.reference_start, "+", "MN908947.3"
                )
                assert lp is not None
                self.assertTrue(
                    record.reference_start >= lp.end
                )  # lp.end is non inclusive

                # Find the right primer
                rp = find_primer_with_lookup(
                    primer_lookup, record.reference_end, "-", "MN908947.3"
                )
                assert rp is not None
                self.assertTrue(
                    record.reference_end <= rp.start  # type: ignore
                )  # record.reference_end is non inclusive

                # Check rg is correct
                if lp.amplicon_number == rp.amplicon_number:
                    rg = str(lp.pool)
                else:
                    rg = "unmatched"
                self.assertEqual(record.get_tag("RG"), rg)

    def test_align_trim_trim_se_norm(self):
        """Tests primers are trimmed correctly"""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-trim_primers_se_norm"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"

            # Create args with test-specific values
            args = create_args(output=output_sam.absolute(), normalise=200)

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)
            amplicons = create_amplicons(scheme.bedlines)
            pools = set([bl.pool for bl in scheme.bedlines])

            ref_lengths = [("MN908947.3", 29903)]
            primer_lookup = create_primer_lookup(ref_lengths, amplicons, 35)

            # Check the out sam is as expected
            for record in pysam.AlignmentFile(str(output_sam), "r"):
                # Find the left primer
                lp = find_primer_with_lookup(
                    primer_lookup, record.reference_start, "+", "MN908947.3"
                )
                assert lp is not None
                self.assertTrue(
                    record.reference_start >= lp.end
                )  # lp.end is non inclusive

                # Find the right primer
                rp = find_primer_with_lookup(
                    primer_lookup, record.reference_end, "-", "MN908947.3"
                )
                assert rp is not None
                self.assertTrue(
                    record.reference_end <= rp.start  # type: ignore
                )  # record.reference_end is non inclusive

                # Check rg is correct
                if lp.amplicon_number == rp.amplicon_number:
                    rg = str(lp.pool)
                else:
                    rg = "unmatched"
                self.assertEqual(record.get_tag("RG"), rg)

    def test_align_trim_write_reports(self):
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-write_report"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"
            report = tempdir_path / "report.tsv"
            amp_depths = tempdir_path / "amp_depths.tsv"

            # Create args with report enabled
            args = create_args(
                output=output_sam.absolute(),
                report=report,
                amp_depth_report=amp_depths,
            )

            # Run
            go(args)
            self.assertTrue(pathlib.Path.exists(report))
            self.assertTrue(pathlib.Path.exists(amp_depths))

            with open(amp_depths, "r") as f:
                next(f)  # Skip header
                for line in f:
                    chrom, amplicon, mean_depth = line.strip().split("\t")
                    self.assertEqual(chrom, "MN908947.3")
                    self.assertIn(
                        amplicon,
                        [
                            str(a.amplicon_number)
                            for a in create_amplicons(
                                Scheme.from_file(args.bedfile).bedlines
                            )
                        ],
                    )
                    self.assertTrue(float(mean_depth) > 0)

    def test_align_trim_require_full_length(self):
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-require_full_length"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"

            # Create args with report enabled
            args = create_args(
                output=output_sam.absolute(),
                require_full_length=True,
                allow_incorrect_pairs=False,
            )

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)
            amplicons = create_amplicons(scheme.bedlines)
            pools = set([bl.pool for bl in scheme.bedlines])

            ref_lengths = [("MN908947.3", 29903)]
            primer_lookup = create_primer_lookup(ref_lengths, amplicons, 35)

            # Check the out sam is as expected
            for record in pysam.AlignmentFile(str(output_sam), "r"):
                # Find the left primer
                lp = find_primer_with_lookup(
                    primer_lookup, record.reference_start, "+", "MN908947.3"
                )
                assert lp is not None
                print(record)
                self.assertLessEqual(
                    record.reference_start,
                    lp.end + 1,
                    "reference_start !<= lp.end",
                )  # lp.end is non inclusive

                # Find the right primer
                rp = find_primer_with_lookup(
                    primer_lookup, record.reference_end, "-", "MN908947.3"
                )
                assert rp is not None
                self.assertGreaterEqual(
                    record.reference_end,
                    rp.start,
                    "reference_end !>= rp.start",
                )  # record.reference_end is non inclusive

    def test_align_trim_se_no_trim(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-se_no_trim") as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"

            # Create args with report enabled
            args = create_args(
                output=output_sam.absolute(),
                require_full_length=True,
                allow_incorrect_pairs=False,
                no_trim_primers=True,
            )

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)
            amplicons = create_amplicons(scheme.bedlines)
            pools = set([bl.pool for bl in scheme.bedlines])

            ref_lengths = [("MN908947.3", 29903)]
            primer_lookup = create_primer_lookup(ref_lengths, amplicons, 35)

            # Check the out sam is as expected
            for record in pysam.AlignmentFile(str(output_sam), "r"):
                # Find the left primer
                lp = find_primer_with_lookup(
                    primer_lookup, record.reference_start, "+", "MN908947.3"
                )
                assert lp is not None
                self.assertLess(
                    record.reference_start,
                    lp.end,
                    f"{record.query_name} - reference_start !< lp.end",
                )  # lp.end is non inclusive

                # Find the right primer
                rp = find_primer_with_lookup(
                    primer_lookup, record.reference_end, "-", "MN908947.3"
                )
                assert rp is not None
                self.assertGreaterEqual(
                    record.reference_end,
                    rp.start,
                    f"{record.query_name} - reference_end !>= rp.start",
                )  # record.reference_end is non inclusive

    def test_align_trim_paired_norm(self):
        with tempfile.TemporaryDirectory(dir="tests", suffix="-paired") as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_bam = tempdir_path / "output.bam"
            amp_depths = tempdir_path / "amp_depths.tsv"
            report = tempdir_path / "report.tsv"

            # Create args with report enabled
            args = create_args(
                output=output_bam.absolute(),
                allow_incorrect_pairs=False,
                amp_depth_report=amp_depths,
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
                normalise=200,
                report=report,
            )

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)

            zero_depth_amps = [15, 29]

            with open(amp_depths, "r") as f:
                next(f)  # Skip header
                for line in f:
                    chrom, amplicon, mean_depth = line.strip().split("\t")
                    self.assertEqual(chrom, "MN908947.3")
                    self.assertIn(
                        amplicon,
                        [
                            str(a.amplicon_number)
                            for a in create_amplicons(
                                Scheme.from_file(args.bedfile).bedlines
                            )
                        ],
                    )
                    if int(amplicon) in zero_depth_amps:
                        self.assertEqual(float(mean_depth), 0.0)
                    else:
                        self.assertTrue(float(mean_depth) > 0)

    def test_align_trim_paired_no_norm(self):
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-paired-no-norm"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_bam = tempdir_path / "output.bam"
            amp_depths = tempdir_path / "amp_depths.tsv"

            # Create args with report enabled
            args = create_args(
                output=output_bam.absolute(),
                allow_incorrect_pairs=False,
                amp_depth_report=amp_depths,
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
            )

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)

            zero_depth_amps = [15, 29]

            with open(amp_depths, "r") as f:
                next(f)  # Skip header
                for line in f:
                    chrom, amplicon, mean_depth = line.strip().split("\t")
                    self.assertEqual(chrom, "MN908947.3")
                    self.assertIn(
                        amplicon,
                        [
                            str(a.amplicon_number)
                            for a in create_amplicons(
                                Scheme.from_file(args.bedfile).bedlines
                            )
                        ],
                    )
                    if int(amplicon) in zero_depth_amps:
                        self.assertEqual(float(mean_depth), 0.0)
                    else:
                        self.assertTrue(float(mean_depth) > 0)

    def test_align_trim_paired_full_length(self):
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-paired-full-length"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_bam = tempdir_path / "output.bam"
            amp_depths = tempdir_path / "amp_depths.tsv"

            # Create args with report enabled
            args = create_args(
                output=output_bam.absolute(),
                allow_incorrect_pairs=False,
                amp_depth_report=amp_depths,
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
                require_full_length=True,
            )

            # Run
            go(args)

            # Read in scheme, create look ups
            scheme = Scheme.from_file(args.bedfile)
            scheme.bedlines = merge_primers(scheme.bedlines)
            amplicons = create_amplicons(scheme.bedlines)
            pools = set([bl.pool for bl in scheme.bedlines])

            ref_lengths = [("MN908947.3", 29903)]
            primer_lookup = create_primer_lookup(ref_lengths, amplicons, 35)

            # Check the out sam is as expected
            for segment1, segment2 in read_pair_generator(
                pysam.AlignmentFile(str(output_bam), "r")
            ):
                if segment1.reference_start < segment2.reference_start:
                    # Find the primers for the segments
                    lp = find_primer_with_lookup(
                        primer_lookup, segment1.reference_start, "+", "MN908947.3"
                    )
                    rp = find_primer_with_lookup(
                        primer_lookup, segment2.reference_end, "-", "MN908947.3"
                    )
                    self.assertLessEqual(
                        segment1.reference_start,
                        lp.end,
                        "forward segment reference_start !<= lp.end",
                    )
                    self.assertGreaterEqual(
                        segment2.reference_end,
                        rp.start,
                        "reverse segment reference_end !>= rp.start",
                    )

                else:
                    # Find the primers for the segments
                    lp = find_primer_with_lookup(
                        primer_lookup, segment2.reference_start, "+", "MN908947.3"
                    )
                    rp = find_primer_with_lookup(
                        primer_lookup, segment1.reference_end, "-", "MN908947.3"
                    )
                    self.assertLessEqual(
                        segment2.reference_start,
                        lp.end,
                        f"{segment2.query_name} - forward segment reference_start !<= lp.end",
                    )
                    self.assertGreaterEqual(
                        segment1.reference_end,
                        rp.start,
                        f"{segment1.query_name} - reverse segment reference_end !>= rp.start",
                    )


if __name__ == "__main__":
    unittest.main()

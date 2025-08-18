import pathlib
import unittest
import argparse
from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import merge_primers
from primalbedtools.amplicons import create_amplicons
import tempfile
import pysam
from align_trim.main import go, create_primer_lookup, find_primer_with_lookup

BED_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/v5.3.2.primer.bed"
BAM_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/sars-cov-2_v5.3.2.bam"


def create_args(**kwargs):
    """Create a fake args object with default values matching main.py argument parser"""
    defaults = {
        "bedfile": str(BED_PATH_V5_3_2.absolute()),
        "bamfile": str(BAM_PATH_V5_3_2.absolute()),
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
    def test_align_trim_trim_primers(self):
        """Tests primers are trimmed correctly"""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-trim_primers"
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
            primer_lookup = create_primer_lookup(ref_lengths, pools, amplicons, 35)

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
                    self.assertTrue(float(mean_depth) >= 0)

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
            primer_lookup = create_primer_lookup(ref_lengths, pools, amplicons, 35)

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


if __name__ == "__main__":
    unittest.main()

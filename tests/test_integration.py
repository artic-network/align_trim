import pathlib
import unittest
import argparse
from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import merge_primers
from primalbedtools.amplicons import create_amplicons
import tempfile
import pysam
from aligntrim.main import go, create_primer_lookup, find_primer_with_lookup

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
        "trim_primers": False,
        "paired": False,
        "no_read_groups": False,
        "verbose": False,
        "remove_incorrect_pairs": False,
        "require_full_length": False,
        "output_sam": None,
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


class TestIntegration(unittest.TestCase):
    def test_aligntrim_trim_primers(self):
        """Tests primers are trimmed correctly"""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-trim_primers"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"

            # Create args with test-specific values
            args = create_args(trim_primers=True, output_sam=str(output_sam.absolute()))

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

    def test_aligntrim_write_report(self):
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-write_report"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"
            report = tempdir_path / "report.csv"

            # Create args with report enabled
            args = create_args(output_sam=str(output_sam.absolute()), report=report)

            # Run
            go(args)
            self.assertTrue(pathlib.Path.exists(report))

    def test_aligntrim_require_full_length(self):
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-require_full_length"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            output_sam = tempdir_path / "output.sam"

            # Create args with report enabled
            args = create_args(
                output_sam=str(output_sam.absolute()),
                require_full_length=True,
                remove_incorrect_pairs=True,
                trim_primers=True,
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
                self.assertEqual(
                    record.reference_start,
                    lp.end + 1 + 1,
                    "reference_start != left primer end",
                )  # lp.end is non inclusive

                # Find the right primer
                rp = find_primer_with_lookup(
                    primer_lookup, record.reference_end, "-", "MN908947.3"
                )
                assert rp is not None
                self.assertEqual(
                    record.reference_end,
                    rp.start - 1,  # type: ignore
                    "reference_end != right primer start",
                )  # record.reference_end is non inclusive


if __name__ == "__main__":
    unittest.main()

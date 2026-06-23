import argparse
import pathlib
import tempfile
import unittest

import pysam

from align_trim.main import go

BED_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/v5.3.2.primer.bed"
BAM_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/sars-cov-2_v5.3.2.bam"
BED_PATH_V3_0_0 = pathlib.Path(__file__).parent / "test_data/v3.0.0.primer.bed"
BAM_PATH_PAIRED_V3_0_0 = (
    pathlib.Path(__file__).parent / "test_data/sars-cov-2_v3.0.0_paired.bam"
)


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
        "genome_coverage_report": None,
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def read_amp_depths(path):
    """Parse an amp_depth_report TSV into {amplicon: mean_depth}."""
    depths = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            _, amplicon, mean_depth = line.strip().split("\t")
            depths[amplicon] = float(mean_depth)
    return depths


def count_reads(sam_path):
    """Count records in a SAM/BAM file."""
    return sum(1 for _ in pysam.AlignmentFile(str(sam_path), "r"))


class TestNormaliseSingleEnd(unittest.TestCase):
    def test_depth_bounded_by_target(self):
        """Mean amplicon depth after SE normalisation should not grossly exceed the target."""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-se_norm_depth_bounded"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            normalise_target = 200

            args = create_args(
                output=(tempdir_path / "output.sam").absolute(),
                normalise=normalise_target,
                amp_depth_report=(tempdir_path / "amp_depths.tsv").absolute(),
            )
            go(args)

            for amplicon, mean_depth in read_amp_depths(
                tempdir_path / "amp_depths.tsv"
            ).items():
                self.assertLessEqual(
                    mean_depth,
                    normalise_target * 1.5,
                    f"Amplicon {amplicon} mean depth {mean_depth:.1f} exceeds 1.5× target {normalise_target}",
                )

    def test_reduces_depth_vs_unnormalized(self):
        """SE normalisation must reduce depth for at least one amplicon vs. unnormalized."""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-se_norm_reduces_depth"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            args_no_norm = create_args(
                output=(tempdir_path / "no_norm.sam").absolute(),
                normalise=0,
                amp_depth_report=(tempdir_path / "depths_no_norm.tsv").absolute(),
            )
            go(args_no_norm)

            # Use a low target (10) to ensure reduction fires even on small test BAMs
            args_norm = create_args(
                output=(tempdir_path / "norm.sam").absolute(),
                normalise=10,
                amp_depth_report=(tempdir_path / "depths_norm.tsv").absolute(),
            )
            go(args_norm)

            no_norm = read_amp_depths(tempdir_path / "depths_no_norm.tsv")
            norm = read_amp_depths(tempdir_path / "depths_norm.tsv")

            reduction_found = any(
                norm[amp] < no_norm[amp]
                for amp in norm
                if no_norm.get(amp, 0) > 0
            )
            self.assertTrue(
                reduction_found,
                "SE normalisation did not reduce depth for any amplicon",
            )

    def test_high_target_keeps_all_reads(self):
        """With a target far above actual coverage, all reads must be accepted."""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-se_norm_high_target"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            args_no_norm = create_args(
                output=(tempdir_path / "no_norm.sam").absolute(),
                normalise=0,
            )
            go(args_no_norm)

            args_high_norm = create_args(
                output=(tempdir_path / "high_norm.sam").absolute(),
                normalise=10000,
            )
            go(args_high_norm)

            count_no_norm = count_reads(tempdir_path / "no_norm.sam")
            count_high_norm = count_reads(tempdir_path / "high_norm.sam")

            self.assertEqual(
                count_no_norm,
                count_high_norm,
                f"High-target SE normalisation dropped reads: "
                f"expected {count_no_norm}, got {count_high_norm}",
            )


class TestNormalisePaired(unittest.TestCase):
    def test_reduces_depth_vs_unnormalized(self):
        """Paired normalisation must reduce depth for at least one amplicon vs. unnormalized."""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-paired_norm_reduces_depth"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            args_no_norm = create_args(
                output=(tempdir_path / "no_norm.bam").absolute(),
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
                normalise=0,
                amp_depth_report=(tempdir_path / "depths_no_norm.tsv").absolute(),
            )
            go(args_no_norm)

            args_norm = create_args(
                output=(tempdir_path / "norm.bam").absolute(),
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
                normalise=200,
                amp_depth_report=(tempdir_path / "depths_norm.tsv").absolute(),
            )
            go(args_norm)

            no_norm = read_amp_depths(tempdir_path / "depths_no_norm.tsv")
            norm = read_amp_depths(tempdir_path / "depths_norm.tsv")

            reduction_found = any(
                norm[amp] < no_norm[amp]
                for amp in norm
                if no_norm.get(amp, 0) > 0
            )
            self.assertTrue(
                reduction_found,
                "Paired normalisation did not reduce depth for any amplicon",
            )

    def test_high_target_keeps_all_reads(self):
        """With a target far above actual coverage, all paired reads must be accepted."""
        with tempfile.TemporaryDirectory(
            dir="tests", suffix="-paired_norm_high_target"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            args_no_norm = create_args(
                output=(tempdir_path / "no_norm.bam").absolute(),
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
                normalise=0,
            )
            go(args_no_norm)

            args_high_norm = create_args(
                output=(tempdir_path / "high_norm.bam").absolute(),
                bedfile=BED_PATH_V3_0_0.absolute(),
                samfile=BAM_PATH_PAIRED_V3_0_0.absolute(),
                normalise=10000,
            )
            go(args_high_norm)

            count_no_norm = count_reads(tempdir_path / "no_norm.bam")
            count_high_norm = count_reads(tempdir_path / "high_norm.bam")

            self.assertEqual(
                count_no_norm,
                count_high_norm,
                f"High-target paired normalisation dropped reads: "
                f"expected {count_no_norm}, got {count_high_norm}",
            )


if __name__ == "__main__":
    unittest.main()

import pathlib
import unittest
from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import BedLine
import pysam

from align_trim.main import find_primer, trim

BED_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/primer.bed"
BED_PATH_V1_0_0 = pathlib.Path(__file__).parent / "test_data/v1.0.0.primer.bed"


bl1 = BedLine("chrom", 0, 10, "primer_1_LEFT_1", 1, "+", "ATCG")
bl2 = BedLine("chrom", 30, 40, "primer_1_RIGHT_1", 1, "-", "ATCG")
bl3 = BedLine("chrom", 10, 20, "primer_2_LEFT_1", 2, "+", "ATCG")
bl4 = BedLine("chrom", 40, 50, "primer_2_RIGHT_1", 2, "-", "ATCG")

dummy_scheme = Scheme(
    None,
    [
        bl1,
        bl2,
        bl3,
        bl4,
    ],
)
# nCov alignment segment (derived from a real nCov read)
seg1 = pysam.AlignedSegment()
seg1.query_name = "0be29940-97ae-440e-b02c-07748edeceec"
seg1.flag = 0
seg1.reference_id = 0
seg1.reference_start = 4294
seg1.mapping_quality = 60
seg1.cigarstring = "40S9M1D8M2D55M1D4M2I12M1I14M1D101M1D53M2D78M1I60M52S"
seg1.query_sequence = "CAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTACAGTGCTTAAAAGTGTAAGTGCCTTTTACATTCTACCATCTATTATCTCTAATGAGAAGCAAGAAATTCTTGAACCTTCATACTTGGAATTTTGCGAGAAATGCTGCACATGCAGAAGAAACACGCAAATTAATGCCTGTCTGTGTGGAAACTAAAGCCATAGTTTCAACTATACAGCGTAAATATAAGGGTATTAAAATACAAGGGGTGTGGTTGATTATGGTGCTAGATTTTACTTTTACACCAGTAAAACAACTGGCGTCACTTATCAACACACTTAACGATCTAAATGAAACTCTTGTTACAATGCCACTTGGCTATGTAACACATGGCTTAGAATTTGGAAGAAGCTGCTCGGTATATGAGATCTCTCAAAGTGCCAGCTACAGTTTCTGTTGCGATTGCTGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACCCATG"
seg1.query_qualities = [30] * 490  # type: ignore

# nCov alignment segment (derived from a real nCov read) - will result in a leading CIGAR deletion during softmasking
seg2 = pysam.AlignedSegment()
seg2.query_name = "15c86d34-a527-4506-9b0c-f62827d01555"
seg2.flag = 0
seg2.reference_id = 0
seg2.reference_start = 4294
seg2.mapping_quality = 60
seg2.cigarstring = (
    "41S9M1D17M1D69M5D1M1D40M2I12M1D41M1D117M1I3M1D4M2D11M2I2M1D18M1D18M1I25M56S"
)
seg2.query_sequence = "CCAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTACAGTGCTTAAAAGTGTAAAAGTGCCTTTACATTCTACCATCTATTATCTCTAATGAGAAGCAAGAAATTCTTGGAACTGTTTCTTGGAATTTGCAGCTTGCACATGCAGAAGAAACACGCAAATTAATGCCTGTCTGTGTGTGGAAACTAAGCCATAGTTTCAACTATACAGCGTAAATATAAGGGTATTAAATACAAGAGGGTGTGGTTGATTATGGTGCTAGATTTTACTTTTACACCAGTAAAACAACTGTAGCGTCACTTATCAACACGCTTAACGATCTAAATGAAACTCTTGTTACAATGCACACTGGCTGTAACACATGAAACTAAATTTGGAAGAAGCTGTCGGTATATGAGATCTCTCCAAAGTGCCAGCTACGGTTTCTGTTAGGTGCTGAAAAGAAAGTTGTCGGTGTCTTTGTGTGAACCTTAGCAATACGTAACC"
seg2.query_qualities = [30] * 490  # type: ignore

# expected softmasked CIGARs
seg1expectedCIGAR = "64S48M1D4M2I12M1I14M1D101M1D53M2D78M1I38M74S"
seg2expectedCIGAR = "67S69M5D1M1D40M2I12M1D41M1D117M1I3M1D4M2D11M2I2M1D18M1D18M1I3M78S"


class TestLegacy(unittest.TestCase):
    """
    Refactor of legacy tests
    """

    def test_find_primer_self(self):
        # test the primer finder on the primers themselves
        for primer in dummy_scheme.bedlines:
            if primer.direction_str == "+":
                result = find_primer(
                    dummy_scheme.bedlines,
                    primer.start,
                    primer.direction_str,
                    primer.chrom,
                )
            else:
                result = find_primer(
                    dummy_scheme.bedlines,
                    primer.end,
                    primer.direction_str,
                    primer.chrom,
                )
            self.assertEqual(result[2].primername, primer.primername)  # type: ignore

    def test_find_primer_other_refs(self):
        # Test other positions
        result = find_primer(
            dummy_scheme.bedlines,
            8,
            "+",
            dummy_scheme.bedlines[0].chrom,
        )
        self.assertEqual(result[2].primername, "primer_2_LEFT_1")  # type: ignore

        result = find_primer(
            dummy_scheme.bedlines,
            25,
            "-",
            dummy_scheme.bedlines[0].chrom,
        )
        self.assertEqual(result[2].primername, "primer_2_LEFT_1")  # type: ignore

    def test_trim_seq(self):
        """test for the trim function"""
        scheme = Scheme.from_file(str(BED_PATH_V1_0_0.absolute()))
        chrom = scheme.bedlines[0].chrom
        data = [(seg1, seg1expectedCIGAR), (seg2, seg2expectedCIGAR)]

        def testRunner(seg, expectedCIGAR):
            # get the nearest primers to the alignment segment
            p1 = find_primer(scheme.bedlines, seg.reference_start, "+", chrom)
            p2 = find_primer(scheme.bedlines, seg.reference_end, "-", chrom)

            # ensure primers were found
            self.assertIsNotNone(
                p1, f"No forward primer found for segment {seg.query_name}"
            )
            self.assertIsNotNone(
                p2, f"No reverse primer found for segment {seg.query_name}"
            )

            # get the primer positions
            p1_position = p1[2].end  # type: ignore
            p2_position = p2[2].start  # type: ignore

            # this segment should need forward and reverse softmasking
            self.assertLess(
                seg.reference_start,
                p1_position,
                f"missed a forward soft masking opportunity (read: {seg.query_name})",
            )
            self.assertGreater(
                seg.reference_end,
                p2_position,
                f"missed a reverse soft masking opportunity (read: {seg.query_name})",
            )

            # before masking, get the query_alignment_length and the CIGAR to use for testing later
            originalCigar = seg.cigarstring
            originalQueryAlnLength = seg.query_alignment_length

            # trim the forward primer
            try:
                trim(seg, p1_position, False, False)
            except Exception as e:
                raise Exception(
                    f"problem soft masking left primer in {seg.query_name} (error: {e})"
                )

            # check the CIGAR and query alignment length is updated
            self.assertNotEqual(
                seg.cigarstring,
                originalCigar,
                f"cigar was not updated with a softmask (read: {seg.query_name})",
            )
            self.assertNotEqual(
                seg.query_alignment_length,
                originalQueryAlnLength,
                f"query alignment was not updated after softmask (read: {seg.query_name})",
            )

            # trim the reverse primer
            try:
                trim(seg, p2_position, True, False)
            except Exception as e:
                raise Exception(
                    f"problem soft masking right primer in {seg.query_name} (error: {e})"
                )

            # check the CIGAR and query alignment length is updated
            self.assertNotEqual(
                seg.cigarstring,
                originalCigar,
                f"cigar was not updated with a softmask (read: {seg.query_name})",
            )
            self.assertNotEqual(
                seg.query_alignment_length,
                originalQueryAlnLength,
                f"query alignment was not updated after softmask (read: {seg.query_name})",
            )

            # check we have the right CIGAR
            self.assertEqual(
                seg.cigarstring,
                expectedCIGAR,
                f"cigar does not match expected cigar string (read: {seg.query_name})",
            )

            # check the query alignment now matches the expected primer product
            self.assertGreaterEqual(
                seg.reference_start,
                p1_position,
                f"left primer not masked correctly (read: {seg.query_name})",
            )
            self.assertLessEqual(
                seg.reference_end,
                p2_position,
                f"right primer not masked correctly (read: {seg.query_name})",
            )

        # run the test with the two alignment segments
        for seg, expectedCIGAR in data:
            with self.subTest(seg=seg.query_name):
                testRunner(seg, expectedCIGAR)


if __name__ == "__main__":
    unittest.main()

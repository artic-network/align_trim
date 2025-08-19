import pathlib
import unittest
from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import merge_primers
from primalbedtools.amplicons import create_amplicons

from align_trim.main import (
    create_primer_lookup,
    find_primer_with_lookup,
)

BED_PATH_V5_3_2 = pathlib.Path(__file__).parent / "test_data/primer.bed"
BED_PATH_V1_0_0 = pathlib.Path(__file__).parent / "test_data/v1.0.0.primer.bed"


class TestCreatePrimerLookup(unittest.TestCase):
    """
    Tests for the create_primer_lookup function
    """

    def setUp(self):
        self.scheme = Scheme.from_file(str(BED_PATH_V5_3_2.absolute()))
        self.pools = {bl.pool for bl in self.scheme.bedlines}
        self.scheme.bedlines = merge_primers(self.scheme.bedlines)
        self.amplicons = create_amplicons(self.scheme.bedlines)

        # Create ref len
        self.ref_len = [("MN908947.3", 29903)]
        return super().setUp()

    def test_create_primer_lookup(self):
        """
        Test that the primer lookup is created correctly
        """
        # Create the primer lookup with 35 padding
        padding = 0
        primer_lookup = create_primer_lookup(
            ref_len_tuple=self.ref_len,
            amplicons=self.amplicons,
            pools=self.pools,
            padding=padding,
        )
        # Check that the primer lookup is a dictionary
        self.assertIsInstance(primer_lookup, dict)
        # Check that the primer lookup contains the expected keys
        self.assertEqual(
            set(x[0] for x in self.ref_len),
            set(primer_lookup.keys()),
        )

        # Check the size of the primer lookup
        self.assertEqual(
            primer_lookup["MN908947.3"].shape,
            (max(self.pools), 29903 + 1),  # +1 for half open interval
        )

        # Check each amplicon is present in the lookup and in correct pool
        for amplicon in self.amplicons:
            # Check amplicon spans correct region, and correct pool
            self.assertEqual(
                set(
                    primer_lookup[amplicon.chrom][
                        amplicon.ipool, amplicon.amplicon_start : amplicon.amplicon_end
                    ]
                ),
                set([amplicon]),
            )
            # Check padding is aligned correctly
            self.assertIsNone(
                primer_lookup[amplicon.chrom][
                    amplicon.ipool, amplicon.amplicon_start - padding - 1
                ]
            )
            self.assertIsNone(
                primer_lookup[amplicon.chrom][
                    amplicon.ipool, amplicon.amplicon_end + padding  # -1 +1
                ]
            )

    def test_create_primer_lookup_padding(self):
        """
        Test that the primer lookup is created correctly
        """
        # Create the primer lookup with 35 padding
        padding = 10
        primer_lookup = create_primer_lookup(
            ref_len_tuple=self.ref_len,
            amplicons=self.amplicons,
            pools=self.pools,
            padding=padding,
        )
        # Check that the primer lookup is a dictionary
        self.assertIsInstance(primer_lookup, dict)
        # Check that the primer lookup contains the expected keys
        self.assertEqual(
            set(x[0] for x in self.ref_len),
            set(primer_lookup.keys()),
        )

        # Check the size of the primer lookup
        self.assertEqual(
            primer_lookup["MN908947.3"].shape,
            (max(self.pools), 29903 + 1),  # +1 for half open interval
        )

        # Check each amplicon is present in the lookup and in correct pool
        for amplicon in self.amplicons:
            # Check amplicon spans correct region, and correct pool
            self.assertEqual(
                set(
                    primer_lookup[amplicon.chrom][
                        amplicon.ipool, amplicon.amplicon_start : amplicon.amplicon_end
                    ]
                ),
                set([amplicon]),
            )
            # Check padding is aligned correctly
            self.assertIsNone(
                primer_lookup[amplicon.chrom][
                    amplicon.ipool, amplicon.amplicon_start - padding - 1
                ]
            )
            self.assertIsNone(
                primer_lookup[amplicon.chrom][
                    amplicon.ipool, amplicon.amplicon_end + padding  # -1 +1
                ]
            )


class TestFindPrimerWithLookup(unittest.TestCase):
    def setUp(self):
        self.scheme = Scheme.from_file(str(BED_PATH_V5_3_2.absolute()))
        self.pools = {bl.pool for bl in self.scheme.bedlines}
        self.scheme.bedlines = merge_primers(self.scheme.bedlines)
        self.amplicons = create_amplicons(self.scheme.bedlines)

        # Create ref len and lookup
        self.ref_len = [("MN908947.3", 29903)]
        self.primer_lookup = create_primer_lookup(
            ref_len_tuple=self.ref_len,
            amplicons=self.amplicons,
            pools=self.pools,
            padding=0,
        )

        return super().setUp()

    def test_find_primer(self):
        # For each position in each amplicons ensure the correct primer is found
        for amplicon in self.amplicons:
            for fpos in range(amplicon.amplicon_start, amplicon.amplicon_start + 50):
                # Find left
                lp = find_primer_with_lookup(
                    self.primer_lookup, fpos, "+", amplicon.chrom
                )
                self.assertEqual(lp, amplicon.left[0])
            for rpos in range(amplicon.amplicon_end - 50, amplicon.amplicon_end):
                # Find Right
                rp = find_primer_with_lookup(
                    self.primer_lookup, rpos, "-", amplicon.chrom
                )
                self.assertEqual(rp, amplicon.right[0])


class TestTrim(unittest.TestCase):
    def setUp(self) -> None:
        return super().setUp()

    def test_trim_basic(self):
        pass


if __name__ == "__main__":
    unittest.main()

import unittest
import time

from ShortSeq import ShortSeq, ShortSeq64, ShortSeq128, ShortSeqVar
from .util import rand_sequence

resources = "./testdata"


class ShortSeqFixedWidthTests(unittest.TestCase):
    """These tests address the fixed-width ShortSeq variants (ShortSeq64 and ShortSeq128)"""

    """Can ShortSeqs represent zero-length sequences? Are they singleton?"""

    def test_empty_seq(self):
        seq_u = ShortSeq.from_str("")
        seq_b = ShortSeq.from_bytes(b"")

        self.assertEqual(seq_b, seq_u)
        self.assertEqual(id(seq_b), id(seq_u))
        self.assertEqual(str(seq_b), "")
        self.assertEqual(str(seq_u), "")

    """Can ShortSeqs encode and decode all valid bases from str object inputs?"""

    def test_single_base_str(self):
        bases = [ShortSeq.from_str(b) for b in "ATGC"]
        self.assertTrue(all(type(b) is ShortSeq64 for b in bases))
        self.assertListEqual([str(b) for b in bases], list("ATGC"))

    """Can ShortSeqs encode and decode all valid bases from bytes object inputs?"""

    def test_single_base_bytes(self):
        to_bytes = lambda x: x.encode()
        bases = [ShortSeq.from_bytes(to_bytes(b)) for b in "ATGC"]
        self.assertTrue(all(type(b) is ShortSeq64 for b in bases))
        self.assertListEqual([str(b) for b in bases], list("ATGC"))

    """Does ShortSeq correctly transition to a larger representative object
    when sequence length crosses the 32 base threshold?"""

    def test_correct_subtype_for_length(self):
        seq_32 = ShortSeq.from_str("A" * 32)
        seq_33 = ShortSeq.from_str("A" * 33)

        self.assertIsInstance(seq_32, ShortSeq64)
        self.assertIsInstance(seq_33, ShortSeq128)

    """Is maximum sequence length correctly enforced?"""

    def test_max_length_exceeded(self):
        max_seq = "ATGC" * 256  # 1024 bases, the maximum allowed
        exc_seq = max_seq + "A"
        no_problem = ShortSeq.from_str(max_seq)
        self.assertEqual(str(no_problem), max_seq)

        with self.assertRaisesRegex(Exception, r"(.*)longer than 1024 bases(.*)"):
            ShortSeq.from_str(exc_seq)

    """Are incompatible sequence characters rejected?"""

    def test_incompatible_seq_chars(self):
        problems_64 =  ["N", "*"]
        problems_128 = [c * 33 for c in problems_64]

        for prob in problems_64 + problems_128:
            with self.assertRaisesRegex(Exception, "Unsupported base character"):
                ShortSeq.from_str(prob)


class ShortSeqVarTests(unittest.TestCase):
    """These tests address the variable length ShortSeq variant (ShortSeqVar)"""

    """Does ShortSeq correctly transition to using ShortSeqVar objects
    when sequence length crosses the 64 base threshold?"""

    def test_min_length(self):
        sample_len = 65
        n_samples = 3

        for _ in range(n_samples):
            sample = rand_sequence(sample_len, no_range=True)
            sq = ShortSeq.from_str(sample)

            self.assertIsInstance(sq, ShortSeqVar)
            self.assertEqual(len(sq), len(sample))
            self.assertEqual(str(sq), sample)

    """Checks that randomly generated sequences encode and decode correctly
    for the entire valid range of lengths."""

    def test_length_range(self):
        # TODO: make importable constants
        min_len, max_len = 65, 1024

        for length in range(min_len, max_len):
            sample = rand_sequence(length, no_range=True)
            sq = ShortSeq.from_str(sample)

            self.assertIsInstance(sq, ShortSeqVar)
            self.assertEqual(len(sq), len(sample))
            self.assertEqual(str(sq), sample)


if __name__ == '__main__':
    unittest.main()

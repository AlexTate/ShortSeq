import unittest

from ShortSeq import ShortSeq, ShortSeq64, ShortSeq128, ShortSeqVar
from ShortSeq import MIN_VAR_NT, MAX_VAR_NT, MIN_64_NT, MAX_64_NT, MIN_128_NT, MAX_128_NT
from ShortSeq.tests.util import rand_sequence, print_var_seq_pext_chunks

resources = "./testdata"


class ShortSeqFixedWidthTests(unittest.TestCase):
    """These tests address the fixed-width ShortSeq variants (ShortSeq64 and ShortSeq128)"""

    """Can ShortSeqs represent zero-length sequences? Are they singleton?"""

    def test_empty_seq(self):
        seq_u = ShortSeq.from_str("")
        seq_b = ShortSeq.from_bytes(b"")

        self.assertEqual(seq_b, seq_u)          # ShortSeq-ShortSeq equality
        self.assertEqual(id(seq_b), id(seq_u))  # singleton
        self.assertEqual(str(seq_b), "")        # decoded string
        self.assertEqual(str(seq_u), "")        # decoded string
        self.assertEqual(seq_b, "")             # __eq__ with bytes argument
        self.assertEqual(seq_u, "")             # __eq__ with str argument

    """Can ShortSeqs encode and decode all valid bases from str object inputs?"""

    def test_single_base_str(self):
        bases = [ShortSeq.from_str(b) for b in "ATGC"]

        self.assertListEqual(bases, list("ATGC"))                    # __eq__ with str argument
        self.assertListEqual([str(b) for b in bases], list("ATGC"))  # decoded string
        self.assertTrue(all(type(b) is ShortSeq64 for b in bases))   # appropriate type

    """Can ShortSeqs encode and decode all valid bases from bytes object inputs?"""

    def test_single_base_bytes(self):
        to_bytes = lambda x: x.encode()
        bases = [ShortSeq.from_bytes(to_bytes(b)) for b in "ATGC"]

        self.assertListEqual(bases, list("ATGC"))                    # __eq__ with str argument
        self.assertListEqual([str(b) for b in bases], list("ATGC"))  # decoded string
        self.assertTrue(all(type(b) is ShortSeq64 for b in bases))   # appropriate type

    """Does ShortSeq correctly transition to a larger representative object
    when sequence length crosses the 32 base threshold?"""

    def test_correct_subtype_for_length(self):
        seq_32 = ShortSeq.from_str("A" * MAX_64_NT)
        seq_33 = ShortSeq.from_str("A" * (MAX_64_NT + 1))

        self.assertIsInstance(seq_32, ShortSeq64)
        self.assertIsInstance(seq_33, ShortSeq128)

    """Are incompatible sequence characters rejected?"""

    def test_incompatible_seq_chars(self):
        problems_64 =  ["N", "*"]
        problems_128 = [c * 33 for c in problems_64]

        for prob in problems_64 + problems_128:
            with self.assertRaisesRegex(Exception, "Unsupported base character"):
                ShortSeq.from_str(prob)

    """Checks that randomly generated sequences encode and decode correctly
        for the entire valid range of lengths."""

    def test_length_range(self):
        # ShortSeq64
        length = None
        try:
            for length in range(MIN_64_NT, MAX_64_NT):
                sample = rand_sequence(length, no_range=True)
                sq = ShortSeq.from_str(sample)

                self.assertIsInstance(sq, ShortSeq64)
                self.assertEqual(len(sq), len(sample))
                self.assertEqual(str(sq), sample)
        except Exception as e:
            print(f"Failed at length {length} (ShortSeq64)")
            raise e

        # ShortSeq128
        length = None
        try:
            for length in range(MIN_128_NT, MAX_128_NT):
                sample = rand_sequence(length, no_range=True)
                sq = ShortSeq.from_str(sample)

                self.assertIsInstance(sq, ShortSeq128)
                self.assertEqual(len(sq), len(sample))
                self.assertEqual(str(sq), sample)
        except Exception as e:
            print(f"Failed at length {length} (ShortSeq128)")
            raise e

    """Can fixed width ShortSeqs be indexed like strings?"""

    def test_subscript(self):
        #ShortSeq64
        sample64, length, i, sq = None, None, None, None
        try:
            for length in range(MIN_64_NT, MAX_64_NT):
                sample64 = rand_sequence(length, no_range=True)
                sq = ShortSeq.from_str(sample64)
                for i in range(len(sample64)):
                    self.assertEqual(sq[i], sample64[i])
                    self.assertEqual(sq[-i], sample64[-i])
        except Exception as e:
            print(f"Failed at length {length} index {i} (ShortSeq64)")
            raise e

        for oob in [len(sample64) + 1, -len(sample64) - 1]:
            with self.assertRaises(IndexError):
                _ = sq[oob]

        # ShortSeq128
        sample128, length, i, sq = None, None, None, None
        try:
            for length in range(MIN_128_NT, MAX_128_NT):
                sample128 = rand_sequence(length, no_range=True)
                sq = ShortSeq.from_str(sample128)
                for i in range(len(sample128)):
                    self.assertEqual(sq[i], sample128[i])
                    self.assertEqual(sq[-i], sample128[-i])
        except Exception as e:
            print(f"Failed at length {length} index {i} (ShortSeq128)")
            raise e

        for oob in [len(sample128) + 1, -len(sample128) - 1]:
            with self.assertRaises(IndexError):
                _ = sq[oob]


    """Can fixed width ShortSeqs be sliced like strings?"""

    def test_slice(self):
        #ShortSeq64
        sample = rand_sequence(MAX_64_NT, no_range=True)
        sq = ShortSeq.from_str(sample)
        self.assertEqual(sq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(sq[:i], sample[:i])
            self.assertEqual(sq[:-i], sample[:-i])
            self.assertEqual(sq[i:], sample[i:])
            self.assertEqual(sq[-i:], sample[-i:])

        # ShortSeq128
        sample = rand_sequence(MAX_128_NT, no_range=True)
        sq = ShortSeq.from_str(sample)
        self.assertEqual(sq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(sq[:i], sample[:i])
            self.assertEqual(sq[:-i], sample[:-i])
            self.assertEqual(sq[i:], sample[i:])
            self.assertEqual(sq[-i:], sample[-i:])


class ShortSeqVarTests(unittest.TestCase):
    """These tests address the variable length ShortSeq variant (ShortSeqVar)"""

    """Does ShortSeq correctly transition to using ShortSeqVar objects
    when sequence length crosses the 64 base threshold?"""

    def test_min_length(self):
        sample_len = MIN_VAR_NT
        n_samples = 3

        for _ in range(n_samples):
            sample = rand_sequence(sample_len, no_range=True)
            sq = ShortSeq.from_str(sample)

            self.assertIsInstance(sq, ShortSeqVar)
            self.assertEqual(len(sq), len(sample))
            self.assertEqual(str(sq), sample)

    """Is maximum sequence length correctly enforced?"""

    def test_max_length(self):
        max_seq = "ATGC" * 256  # 1024 bases, the maximum allowed
        exc_seq = max_seq + "A"
        no_problem = ShortSeq.from_str(max_seq)
        self.assertEqual(str(no_problem), max_seq)

        with self.assertRaisesRegex(Exception, r"(.*)longer than 1024 bases(.*)"):
            ShortSeq.from_str(exc_seq)

    """Checks that randomly generated sequences encode and decode correctly
    for the entire valid range of lengths."""

    def test_length_range(self):
        length = None
        try:
            for length in range(MIN_VAR_NT, MAX_VAR_NT-1):
                sample = rand_sequence(length, no_range=True)
                sq = ShortSeq.from_str(sample)

                self.assertIsInstance(sq, ShortSeqVar)
                self.assertEqual(len(sq), len(sample))
                self.assertEqual(str(sq), sample)
        except Exception as e:
            print(f"Failed at length {length}")
            raise e

    """Can ShortSeqVars be indexed like strings?"""

    def test_subscript(self):
        length, i = None, None
        try:
            for length in range(MIN_VAR_NT, MAX_VAR_NT-1):
                sample = rand_sequence(length, no_range=True)
                sq = ShortSeq.from_str(sample)
                for i in range(len(sample)):
                    self.assertEqual(sq[i], sample[i])
                    self.assertEqual(sq[-i], sample[-i])
        except Exception as e:
            print(f"Failed at length {length} index {i}")
            raise e

        for oob in [len(sample) + 1, -len(sample) - 1]:
            with self.assertRaises(IndexError):
                _ = sq[oob]

    """Can ShortSeqVars be sliced like strings?"""

    def test_slice(self):
        # Min length
        sample = rand_sequence(MIN_VAR_NT, no_range=True)
        sq = ShortSeq.from_str(sample)
        self.assertEqual(sq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(sq[:i], sample[:i])
            self.assertEqual(sq[:-i], sample[:-i])
            self.assertEqual(sq[i:], sample[i:])
            self.assertEqual(sq[-i:], sample[-i:])

        # Max length
        sample = rand_sequence(MAX_VAR_NT, no_range=True)
        sq = ShortSeq.from_str(sample)
        self.assertEqual(sq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(sq[:i], sample[:i])
            self.assertEqual(sq[:-i], sample[:-i])
            self.assertEqual(sq[i:], sample[i:])
            self.assertEqual(sq[-i:], sample[-i:])


if __name__ == '__main__':
    unittest.main()

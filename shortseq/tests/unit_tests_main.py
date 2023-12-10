import unittest
import sys

from random import randint

import shortseq as sq
from shortseq import ShortSeq64, ShortSeq128, ShortSeqVar
from shortseq import MIN_VAR_NT, MAX_VAR_NT, MIN_64_NT, MAX_64_NT, MIN_128_NT, MAX_128_NT
from shortseq.tests.util import rand_sequence, print_var_seq_pext_chunks

resources = "./testdata"


class ShortSeqFixedWidthTests(unittest.TestCase):
    """These tests address the fixed-width ShortSeq variants (ShortSeq64 and ShortSeq128)"""

    """Can ShortSeqs represent zero-length sequences? Are they singleton?"""

    def test_empty_seq(self):
        seq_u = sq.pack("")
        seq_b = sq.pack(b"")

        self.assertEqual(seq_b, seq_u)          # ShortSeq-ShortSeq equality
        self.assertEqual(id(seq_b), id(seq_u))  # singleton
        self.assertEqual(str(seq_b), "")        # decoded string
        self.assertEqual(str(seq_u), "")        # decoded string
        self.assertEqual(seq_b, "")             # __eq__ with bytes argument
        self.assertEqual(seq_u, "")             # __eq__ with str argument

    """Can ShortSeqs encode and decode all valid bases from str object inputs?"""

    def test_single_base_str(self):
        bases = [sq.from_str(b) for b in "ATGC"]

        self.assertListEqual(bases, list("ATGC"))                    # __eq__ with str argument
        self.assertListEqual([str(b) for b in bases], list("ATGC"))  # decoded string
        self.assertTrue(all(type(b) is ShortSeq64 for b in bases))   # appropriate type

    """Can ShortSeqs encode and decode all valid bases from bytes object inputs?"""

    def test_single_base_bytes(self):
        to_bytes = lambda x: x.encode()
        bases = [sq.from_bytes(to_bytes(b)) for b in "ATGC"]

        self.assertListEqual(bases, list("ATGC"))                    # __eq__ with str argument
        self.assertListEqual([str(b) for b in bases], list("ATGC"))  # decoded string
        self.assertTrue(all(type(b) is ShortSeq64 for b in bases))   # appropriate type

    """Does ShortSeq correctly transition to a larger representative object
    when sequence length crosses the 32 base threshold?"""

    def test_correct_subtype_for_length(self):
        seq_32 = sq.pack("A" * MAX_64_NT)
        seq_33 = sq.pack("A" * (MAX_64_NT + 1))

        self.assertIsInstance(seq_32, ShortSeq64)
        self.assertIsInstance(seq_33, ShortSeq128)

    """Are incompatible sequence characters rejected?"""

    def test_incompatible_seq_chars(self):
        problems_64 =  ["N", "*"]
        problems_128 = [c * 33 for c in problems_64]

        for prob in problems_64 + problems_128:
            with self.assertRaisesRegex(Exception, "Unsupported base character"):
                sq.pack(prob)

    """Are min and max length ShortSeqs the correct size?"""

    def test_size(self):
        # ShortSeq64
        seq_min = sq.pack(rand_sequence(MIN_64_NT))
        seq_max = sq.pack(rand_sequence(MAX_64_NT))

        self.assertEqual(sys.getsizeof(seq_min), 32)
        self.assertEqual(sys.getsizeof(seq_max), 32)

        # ShortSeq128
        seq_min = sq.pack(rand_sequence(MIN_128_NT))
        seq_max = sq.pack(rand_sequence(MAX_128_NT))

        self.assertEqual(sys.getsizeof(seq_min), 48)
        self.assertEqual(sys.getsizeof(seq_max), 48)

    """Checks that randomly generated sequences encode and decode correctly
        for the entire valid range of lengths."""

    def test_length_range(self):
        # ShortSeq64
        length = None
        try:
            for length in range(MIN_64_NT, MAX_64_NT):
                sample = rand_sequence(length)
                seq = sq.pack(sample)

                self.assertIsInstance(seq, ShortSeq64)
                self.assertEqual(len(seq), len(sample))
                self.assertEqual(str(seq), sample)
        except Exception as e:
            print(f"Failed at length {length} (ShortSeq64)")
            raise e

        # ShortSeq128
        length = None
        try:
            for length in range(MIN_128_NT, MAX_128_NT):
                sample = rand_sequence(length)
                seq = sq.pack(sample)

                self.assertIsInstance(seq, ShortSeq128)
                self.assertEqual(len(seq), len(sample))
                self.assertEqual(str(seq), sample)
        except Exception as e:
            print(f"Failed at length {length} (ShortSeq128)")
            raise e

    """Can fixed width ShortSeqs be indexed like strings?"""

    def test_subscript(self):
        #ShortSeq64
        sample64, length, i, seq = None, None, None, None
        try:
            for length in range(MIN_64_NT, MAX_64_NT):
                sample64 = rand_sequence(length)
                seq = sq.pack(sample64)
                for i in range(len(sample64)):
                    self.assertEqual(seq[i],  sample64[i])
                    self.assertEqual(seq[-i], sample64[-i])
        except Exception as e:
            print(f"Failed at length {length} index {i} (ShortSeq64)")
            raise e

        for oob in [len(sample64) + 1, -len(sample64) - 1]:
            with self.assertRaises(IndexError):
                _ = seq[oob]

        # ShortSeq128
        sample128, length, i, seq = None, None, None, None
        try:
            for length in range(MIN_128_NT, MAX_128_NT):
                sample128 = rand_sequence(length)
                seq = sq.pack(sample128)
                for i in range(len(sample128)):
                    self.assertEqual(seq[i],  sample128[i])
                    self.assertEqual(seq[-i], sample128[-i])
        except Exception as e:
            print(f"Failed at length {length} index {i} (ShortSeq128)")
            raise e

        for oob in [len(sample128) + 1, -len(sample128) - 1]:
            with self.assertRaises(IndexError):
                _ = seq[oob]

    """Does the Hamming distance between two ShortSeqs work as expected?"""

    def test_hamming_distance(self):
        def str_ham(a, b): return sum(a_nt != b_nt for a_nt, b_nt in zip(a, b))

        for length in range(0, MAX_128_NT):
            a = rand_sequence(length)
            b = rand_sequence(length)

            self.assertEqual(sq.pack(a) ^ sq.pack(b), str_ham(a, b))

    """Can fixed width ShortSeqs be sliced like strings?"""

    def test_slice(self):
        #ShortSeq64
        sample = rand_sequence(MAX_64_NT)
        seq = sq.pack(sample)
        self.assertEqual(seq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(seq[:i],  sample[:i])
            self.assertEqual(seq[:-i], sample[:-i])
            self.assertEqual(seq[i:],  sample[i:])
            self.assertEqual(seq[-i:], sample[-i:])

        # ShortSeq128
        sample = rand_sequence(MAX_128_NT)
        seq = sq.pack(sample)
        self.assertEqual(seq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(seq[:i],  sample[:i])
            self.assertEqual(seq[:-i], sample[:-i])
            self.assertEqual(seq[i:],  sample[i:])
            self.assertEqual(seq[-i:], sample[-i:])


class ShortSeqVarTests(unittest.TestCase):
    """These tests address the variable length ShortSeq variant (ShortSeqVar)"""

    """Does ShortSeq correctly transition to using ShortSeqVar objects
    when sequence length crosses the 64 base threshold?"""

    def test_min_length(self):
        sample_len = MIN_VAR_NT
        n_samples = 3

        for _ in range(n_samples):
            sample = rand_sequence(sample_len)
            seq = sq.pack(sample)

            self.assertIsInstance(seq, ShortSeqVar)
            self.assertEqual(len(seq), len(sample))
            self.assertEqual(str(seq), sample)

    """Is maximum sequence length correctly enforced?"""

    def test_max_length(self):
        max_seq = "ATGC" * 256  # 1024 bases, the maximum allowed
        exc_seq = max_seq + "A"
        no_problem = sq.pack(max_seq)
        self.assertEqual(str(no_problem), max_seq)

        with self.assertRaisesRegex(Exception, r"(.*)longer than 1024 bases(.*)"):
            sq.pack(exc_seq)

    """Checks that randomly generated sequences encode and decode correctly
    for the entire valid range of lengths."""

    def test_length_range(self):
        length = None
        try:
            for length in range(MIN_VAR_NT, MAX_VAR_NT):
                sample = rand_sequence(length)
                seq = sq.pack(sample)

                self.assertIsInstance(seq, ShortSeqVar)
                self.assertEqual(len(seq), len(sample))
                self.assertEqual(str(seq), sample)
        except Exception as e:
            print(f"Failed at length {length}")
            raise e

    """Can ShortSeqVars be indexed like strings?"""

    def test_subscript(self):
        length, i = None, None
        try:
            for length in range(MIN_VAR_NT, MAX_VAR_NT):
                sample = rand_sequence(length)
                seq = sq.pack(sample)
                for i in range(len(sample)):
                    self.assertEqual(seq[i],  sample[i])
                    self.assertEqual(seq[-i], sample[-i])
        except Exception as e:
            print(f"Failed at length {length} index {i}")
            raise e

        for oob in [len(sample) + 1, -len(sample) - 1]:
            with self.assertRaises(IndexError):
                _ = seq[oob]

    """Can ShortSeqVars be sliced like strings?"""

    def test_slice(self):
        # Min length
        sample = rand_sequence(MIN_VAR_NT)
        seq = sq.pack(sample)
        self.assertEqual(seq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(seq[:i],  sample[:i])
            self.assertEqual(seq[:-i], sample[:-i])
            self.assertEqual(seq[i:],  sample[i:])
            self.assertEqual(seq[-i:], sample[-i:])

        # Max length
        sample = rand_sequence(MAX_VAR_NT)
        seq = sq.pack(sample)
        self.assertEqual(seq[:], sample)
        for i in range(len(sample)):
            self.assertEqual(seq[:i],  sample[:i])
            self.assertEqual(seq[:-i], sample[:-i])
            self.assertEqual(seq[i:],  sample[i:])
            self.assertEqual(seq[-i:], sample[-i:])

    """Just slice the heck out of the darn thing"""

    def test_stochastic_slice(self):
        sample = rand_sequence(MAX_VAR_NT)
        seq = sq.pack(sample)

        for _ in range(10000):
            a = randint(0, MAX_VAR_NT // 2)
            b = randint(a, a + randint(1, MAX_VAR_NT - a))
            self.assertEqual(seq[a:b], sample[a:b])

    """Does the Hamming distance between two ShortSeqs work as expected?"""

    def test_hamming_distance(self):
        def str_ham(a, b): return sum(a_nt != b_nt for a_nt, b_nt in zip(a, b))

        for length in range(MIN_VAR_NT, MAX_VAR_NT):
            a = rand_sequence(length)
            b = rand_sequence(length)

            self.assertEqual(sq.pack(a) ^ sq.pack(b), str_ham(a, b))

    def test_readme(self):
        # Construct from PyUnicode or PyBytes
        seq_str = "ATGC"
        seq_bytes = b"ATGC"
        seq_1 = sq.pack(seq_str)
        seq_2 = sq.pack(seq_bytes)

        # Verify outputs (optional)
        assert seq_1 == seq_2 == seq_str
        assert len(seq_1) == len(seq_2) == len(seq_str)

        seq_3 = sq.pack("TATTAGCGATTGACAGTTGTCCTGTAATAACGCCGGGTAAATTTGCCG")
        seq_4 = sq.pack("TATTACCGATTGACAGTTGTCCTGTAATAACGGCGGGTAAATTTGCTG")  # 5M1X26M1X13M1X1M
        seq_str = str(seq_4)

        # Slice and subscript
        assert seq_4[5:15] == seq_str[5:15]
        assert seq_4[-2] == seq_str[-2]

        # Vectorized hamming distance (differing bases)
        hammd = sum(a!=b for a, b in zip(seq_3, seq_4))
        assert seq_3 ^ seq_4 == hammd == 3

        # Count unique sequences similar to collections.Counter
        from shortseq import ShortSeqCounter
        counts = ShortSeqCounter([seq_bytes] * 10)
        assert counts == {sq.pack("ATGC"): 10}

    """Are min and max length ShortSeqVars the correct size?"""

    def test_size(self):
        seq_min = sq.pack(rand_sequence(MIN_VAR_NT))
        seq_max = sq.pack(rand_sequence(MAX_VAR_NT))

        self.assertEqual(sys.getsizeof(seq_min), 56)
        self.assertEqual(sys.getsizeof(seq_max), 288)



if __name__ == '__main__':
    unittest.main()

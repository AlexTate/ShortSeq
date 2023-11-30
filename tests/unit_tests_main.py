import unittest
import numpy as np
import psutil
import random
import math
import time
import os

from collections import Counter, defaultdict
from random import randint
from pympler.asizeof import asizeof

from ShortSeq import ShortSeqCounter, read_and_count_fastq, ShortSeq, ShortSeq64, ShortSeq128, ShortSeqVar

resources = "./testdata"

# === HELPER FUNCTIONS ================================================================================

def print_var_seq_pext_chunks(seq):
    """Prints the specified sequence with the block and pext boundaries of ShortSeqVar indicated.
    This is useful for visualizing and debugging the ShortSeqVar marshalling algorithm."""

    block_count = math.ceil(len(seq) / 32)
    blocks = [seq[i*32:(i+1)*32] for i in range(block_count)]
    out = []

    for block in blocks:
        np64, rem = divmod(len(block), 8)
        np32, ser = divmod(rem, 4)
        chunks = []

        if np64: chunks.extend(block[i*8:(i+1)*8] for i in range(np64))
        if np32: chunks.append(block[np64*8:np64*8+4])
        if ser: chunks.append(block[-1 * ser:])

        out.append("|".join(chunks))

    print(" -> ".join(out))


def rand_sequence(min_length=None, max_length=None, no_range=False, as_bytes=False):
    """Returns a randomly generated sequence of the specified type, with a length in the specified range"""

    assert (min_length, max_length) != (None, None)

    if no_range:
        max_length = min_length
        min_length = 0
    if min_length is None:
        min_length = max_length
    if max_length is None:
        max_length = min_length

    seq = ''.join(np.random.choice(["A", "C", "T", "G"]) for _ in range(min_length, max_length))
    return seq.encode() if as_bytes else seq


# === TESTS ===========================================================================================

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

    """A primitive benchmark for PyBytes to ShortSeqVar vs PyUnicode runtime."""

    @unittest.skip("Long running benchmark test")
    def test_benchmark(self):
        sample_size = 10000
        min_len, max_len = 65, 1023
        samples = [rand_sequence(min_len, max_len, as_bytes=True) for _ in range(sample_size)]
        print("samples done")

        start = time.time()
        uni = [seq.decode() for seq in samples]
        end = time.time()
        print(f"{end - start:.8f}s unicode")

        start = time.time()
        sq = [ShortSeq.from_bytes(seq) for seq in samples]
        end = time.time()
        print(f"{end - start:.8f}s ShortSeq")


# For now, you must first generate these files using the make_data() helper function.
# These files will be several hundred megabytes in size.
ten_m_wo_repeat = f"{resources}/10m_no_rep.txt"
one_m_w_repeats = f"{resources}/1m.txt"


class Profiling(unittest.TestCase):
    def make_data(self, n, file=f"{resources}/test_data.txt", repeats=False):
        dups = lambda: randint(1, 100) if repeats else lambda: 1
        mapr = ["A", "T", "G", "C"]
        seqs = []

        for i in range(n):
            length = randint(15, 32)
            new_seq = [mapr[randint(0, 3)] for _ in range(length)]
            seqs.extend(["".join(new_seq)] * dups())

        with open(file, 'w') as f:
            f.write('\n'.join(seqs))

        print("Wrote %d sequences to %s" % (len(seqs), file))

    def read_data(self, file="test_data.txt", n=None):
        with open(file, 'rb') as f:
            seqs = f.readlines()

        length = n if n is not None else len(seqs)
        result = [seqs[i][:-1] for i in range(length)]

        return result

    def benchmark_unicode_tp_hash(self, seqs, decode=False):
        start = time.time()
        count = Counter(seqs)
        end = time.time()
        diff = end - start
        print(f"Time to count bytes: {diff:.6f}s")

        if decode:
            start2 = time.time()
            for seq in count.keys():
                seq.decode()
            end2 = time.time()
            diff = end2 - start
            print(f"Time to decode: {end2 - start2:.6f}")
            print(f"Unicode total: {diff:.6f}")

        print(f"{len(count)} unique sequences")
        return count, diff

    def benchmark_ShortSeq(self, seqs, unmarshall=False):
        start = time.time()
        count = ShortSeqCounter(seqs)
        end = time.time()
        diff = end - start
        print(f"Time to count ShortSeqs: {diff:.6f}")

        if unmarshall:
            start2 = time.time()
            for seq in count.keys():
                seq.__str__()
            end2 = time.time()
            diff = end2 - start
            print(f"Time to unmarshall: {end2 - start2:.6f}")
            print(f"ShortSeq total: {diff:.6f}")

        print(f"{len(count)} unique sequences")
        return count, diff

    def print_compression_ratio(self, old_ctr, new_ctr):
        t1, t2 = 0, 0
        for c1, c2 in zip(old_ctr.keys(), new_ctr.keys()):
            t1 += asizeof(c1)
            t2 += asizeof(c2)

        print(f"Compression Ratio: {t2 / t1}")

    def print_time_diff(self, old_t, new_t):
        fastest = min(old_t, new_t)
        slowest = max(old_t, new_t)
        diff = slowest - fastest
        if new_t == fastest:
            direction = "faster"
            pct = (diff / slowest) * 100
        else:
            direction = "slower"
            pct = (diff / fastest) * 100

        print(f"Time difference: {diff:.5f} ({pct:.1f}% {direction})")

    @unittest.skip("Long running benchmark test")
    def test_benchmark(self):
        mem = []
        process = psutil.Process(os.getpid())
        t1 = time.time()
        # counts = read_and_count_fastq(one_m_w_repeats)
        counts = read_and_count_fastq(ten_m_wo_repeat)
        t2 = time.time()
        print(t2 - t1)
        # sys.exit(0)

        seqs = self.read_data(ten_m_wo_repeat)
        # seqs = self.read_data(one_m_w_repeats)
        mem.append(process.memory_info().rss)

        counts1, t1 = self.benchmark_unicode_tp_hash(seqs, decode=True)
        mem.append(process.memory_info().rss)
        print('=' * 35)
        counts2, t2 = self.benchmark_ShortSeq(seqs, unmarshall=True)

        mem.append(process.memory_info().rss)
        print('=' * 35)
        print(f"Size of PyUnicode dictionary: {asizeof(counts1)/1024/1024:.1f}mb")
        print(f"Size of ShortSeqCounter dict: {asizeof(counts2)/1024/1024:.1f}mb")
        print([f"{int(m / 1024 / 1024):d}mb" for m in mem])
        print('=' * 35)

        self.print_compression_ratio(old_ctr=counts1, new_ctr=counts2)
        self.print_time_diff(t1, t2)
        print(f"Count values match: {sorted(counts1.values()) == sorted(counts2.values())}")

if __name__ == '__main__':
    unittest.main()

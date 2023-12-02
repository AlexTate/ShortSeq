import numpy as np
import math


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
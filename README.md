# ShortSeq

ShortSeqs are compact and efficient Python objects that hold short sequences while using up to 73% less memory compared to built-in types. They are prehashed and comparable, they support slicing and indexing, and they easily convert back to their original string form.

| Sequence Length | PyUnicode Size<sup>*</sup> | PyBytes Size<sup>*</sup> |  ShortSeq Size<sup>*</sup> | % Mem. Reduction |
|-----------------|----------------------------|--------------------------|---------------------------:|------------------|
| 0-32 nt         | 56-88 bytes                | 40-72 bytes              |           32 bytes (fixed) | **20-64%**       |
| 33-64 nt        | 88-120 bytes               | 72-104 bytes             |           48 bytes (fixed) | **33-60%**       |
| 65-1024 nt      | 120-1080 bytes             | 104-1064 bytes           |               48-288 bytes | **53-73%**       |

<sup>* Object sizes were measured on Python 3.10 using `asizeof()` from the `pympler` package.</sup>

In the table above, you can see that Python's memory representation of DNA sequences is larger than a C-style `char *` array, which would only need one byte per base. Using Cython we can move some of this memory representation out of Python space and into C space for faster facilities and a more compact bitwise representation.  


### Installation

```shell
mamba install -c bioconda -c conda-forge shortseq
```


### Usage

```python
from ShortSeq import ShortSeq, ShortSeqCounter

# Construct from PyUnicode
seq_str = "ATGC"
seq_1 = ShortSeq.from_str(seq_str)

# Or, construct from PyBytes
seq_bytes = b"ATGC"
seq_2 = ShortSeq.from_bytes(seq_bytes)

# Verify outputs (optional)
assert seq_1 == seq_2
assert seq_str == str(seq_1) == str(seq_2)
assert len(seq_str) == len(seq_1) == len(seq_2)

# Count unique sequences
counts = ShortSeqCounter([seq_bytes] * 10)
assert counts == {ShortSeq.from_str("ATGC"): 10}
```

### CPU Requirements

- Intel Haswell (2014) and newer, or
- AMD Excavator (2015) and newer, or
- Apple M1 and newer

However, AMD processors [prior to Zen 3](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#cite_ref-12) (2020) aren't recommended for 65-1024 nt sequences if runtime performance is a high priority.


### Encoding (Compression)

We represent DNA with four symbols: A, C, T, and G. Generally speaking, when these symbols are represented in computer systems, each symbol takes up one byte or 8 bits of memory because these letters are part of a symbol system that requires 7/8 of those bits for its range. These symbols can instead be represented ordinally as 0, 1, 2, and 3, which only requires 2 bits, allowing us to pack 4 nucleotides into each byte rather than just one.

| Nucleotide | In ASCII     |  Ordinal Value  |
|------------|--------------|:---------------:|
| A          | `0100 00 01` |      `00`       |
| C          | `0100 01 11` |      `01`       |
| T          | `0101 10 00` |      `10`       |
| G          | `0100 11 11` |      `11`       |

This table shows how each nucleotide is represented in ASCII and how it's ordinal value is represented in binary. I should mention that this scheme isn't my work, but rather a well-known technique that's been around for a while.


### Decoding

ShortSeqs are decoded back to their original sequence strings "lazily", i.e. it happens only when you ask and the result isn't cached in the object, so it has to be recomputed with each request. However, ShortSeqs retain the original string's length and can be compared to each other for equality **without** decoding.


### Overhead

The time it takes to convert a list of PyBytes sequences to a list of ShortSeqs is roughly the same as converting the list to PyUnicode sequences. The time it takes to count unique sequences in each of these lists (ShortSeq vs. PyUnicode) is also roughly the same.


### Longer Sequences, Featuring SIMD Encoding

Longer sequences (65 - 1024 bases) are experimentally supported. These sequences are encoded using SIMD instructions which convert 8 bases at a time for higher throughput. However, sequences of this length are not checked for invalid base characters; it is the user's responsibility to ensure that valid DNA strings are used.


### Acknowledgements

Huge thanks to the Montgomery Lab at Colorado State University. This was an experiment of mine that was started in pursuit of faster sequence deduplication and optimization while working on [tinyRNA](https://www.github.com/MontgomeryLab/tinyRNA).

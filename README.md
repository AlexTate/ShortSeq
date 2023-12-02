# ShortSeq

ShortSeqs are compact and efficient Python objects that hold short sequences while using up to 73% less memory compared to built-in types. They have a pre-computed hash value, can be compared for equality, and are easily converted back to the original sequence string.


| Sequence Length | PyUnicode Size<sup>*</sup> | PyBytes Size<sup>*</sup> | ShortSeq Size<sup>*</sup> | % Reduction |
|-----------------|----------------------------|--------------------------|---------------------------|-------------|
| 0-32 nt         | 56-88 bytes                | 40-72 bytes              | 32 bytes (fixed)          | **20-64%**  |
| 33-64 nt        | 88-120 bytes               | 72-104 bytes             | 48 bytes (fixed)          | **33-60%**  |
| 65-1024 nt      | 120-1080 bytes             | 104-1064 bytes           | 48-288 bytes              | **53-73%**  |

<sup>* Object sizes were measured on Python 3.10 using `asizeof()` from the `pympler` package.</sup>

### CPU Requirements

Your processor must support BMI2 instructions. Roughly speaking, this includes:
- Intel Haswell (2014) and newer
- AMD Excavator (2015) and newer
- Apple M1 and newer

However, AMD processors [prior to Zen 3](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#cite_ref-12) (2020) aren't recommended for 65-1024 nt sequences if runtime performance is a high priority.

### Usage
```python
from ShortSeq import ShortSeq, ShortSeqCounter

# Construct from PyUnicode
seq_str = "ATGC"
seq_1 = ShortSeq.from_str(seq_str)

# Construct from PyBytes
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

### Encoding (Compression)

We represent DNA with four symbols: A, C, T, and G. Generally speaking, when these symbols are represented in computer systems, each symbol takes up one byte or 8 bits of memory because these letters are part of a symbol system that requires 7/8 of those bits for its range. These symbols can instead be represented ordinally as 0, 1, 2, and 3, which only requires 2 bits, allowing us to pack 4 nucleotides into each byte rather than just one.

| Nucleotide | In ASCII     |  Ordinal Value  |
|------------|--------------|:---------------:|
| A          | `0100 00 01` |      `00`       |
| C          | `0100 01 11` |      `01`       |
| T          | `0101 10 00` |      `10`       |
| G          | `0100 11 11` |      `11`       |

When Python holds sequences of DNA in memory as PyUnicode or PyBytes objects, their footprint is even larger and more complex than the ASCII representation shown above. Using Cython we can move this memory representation out of Python space and into C space for better efficiency while still fulfilling the duties of a Python object.

### Decoding
ShortSeqs are decoded back to their original sequence strings "lazily", i.e. it happens only when you ask and the result isn't cached in the object, so it has to be recomputed with each request. However, ShortSeqs retain the original string's length and ability to be compared for equality without conversion. The packed integer form of each sequence additionally acts as its pre-computed hash value. These two properties make ShortSeqs ideal for use in datastructures like sets and dictionaries, making them useful for efficiently gathering statistics involving very large volumes of short sequences.

### Overhead
The time it takes to convert a list of PyBytes sequences to a list of ShortSeqs is roughly the same as converting the list to PyUnicode sequences. The time it takes to count unique sequences in each of these lists (ShortSeq vs. PyUnicode) is also roughly the same.

### Longer Sequences, Featuring SIMD Encoding
Longer sequences (65 - 1024 bases) are experimentally supported. These sequences are encoded using SIMD instructions which convert 8 bases at a time for higher throughput. However, sequences of this length are not checked for invalid base characters; it is the user's responsibility to ensure that valid DNA strings are used.
# ShortSeq
Compact representation of short sequences using bit packed integers.

### Compression

We represent DNA with four symbols: A, T, G, and C. Generally speaking, when these symbols are represented in computer systems, each symbol takes up one byte or 8 bits of memory because these letters are part of a symbol system that requires all 8 of those bits for its range. These symbols can instead be represented ordinally as 0, 1, 2, and 3, which only requires 2 bits, so we pack 4 values per byte.

| Nucleotide | In ASCII     |  Ordinal Value  |
|------------|--------------|:---------------:|
| A          | `0100 00 01` |      `00`       |
| C          | `0100 01 11` |      `01`       |
| T          | `0101 10 00` |      `10`       |
| G          | `0100 11 11` |      `11`       |

Since Python is a loosely typed language with a lot of conveniences, both the numerical and string encodings of these symbols are much larger and more complex in their memory representation. 

| Sequence Length | PyUnicode Size<sup>*</sup> | PyBytes Size<sup>*</sup> | ShortSeq Size<sup>*</sup> | % Reduction (vs. PyUnicode) |
|-----------------|----------------------------|--------------------------|---------------------------|-----------------------------|
| 33-64 nt        | 88-120 bytes               | 72-104 bytes             | 48 bytes (fixed)          | 33-60%                      |
| 0-32 nt         | 56-88 bytes                | 40-72 bytes              | 32 bytes (fixed)          | 20-64%                      |

<sup>*</sup> Object sizes were measured on Python 3.10 using `asizeof()` from the `pympler` package.

Using Cython we can move this memory representation out of Python space and into C space for maximum efficiency while still fulfilling the duties of a Python object. A certain amount of space goes to the Python object header and padding to 16 byte multiples for proper memory alignment, so the actual memory savings are less than the theoretical 75% reduction.

### Overhead
The time it takes to convert a list of Python bytes objects (each representing the ASCII encoding of the sequence) to a list of ShortSeq objects is roughly the same. The time it takes to count the number of occurrences of each sequence in each of those lists is also roughly the same. Usually there is 5-15% runtime overhead but the current implementation is sensitive to CPU cache residency, so its range is more variable than that of PyUnicode.

### Decoding
The original ASCII form of each ShortSeq is decoded "lazily", i.e. it happens only when asked and the result isn't cached in the object, so it must be recomputed with each request. However, ShortSeqs still retain the original string's length and their ability to be compared for equality. The packed integer form of each sequence additionally acts as its pre-computed hash value. These two properties make ShortSeqs ideal for use in datastructures like sets and dictionaries, making them useful for efficiently gathering statistics involving very large volumes of short sequences.
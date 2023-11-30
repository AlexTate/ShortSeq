# ShortSeq
Compact representation of short sequences using bit packed integers.

### Compression

We represent DNA with four symbols: A, C, T, and G. Generally speaking, when these symbols are represented in computer systems, each symbol takes up one byte or 8 bits of memory because these letters are part of a symbol system that requires 7/8 of those bits for its range. These symbols can instead be represented ordinally as 0, 1, 2, and 3, which only requires 2 bits, allowing us to represent 4 nucleotides per byte rather than just one.

| Nucleotide | In ASCII     |  Ordinal Value  |
|------------|--------------|:---------------:|
| A          | `0100 00 01` |      `00`       |
| C          | `0100 01 11` |      `01`       |
| T          | `0101 10 00` |      `10`       |
| G          | `0100 11 11` |      `11`       |

Since Python is a loosely typed language with a lot of conveniences, both the numerical and string encodings of these symbols are much larger and more complex in their memory representation. 

| Sequence Length | PyUnicode Size<sup>*</sup> | PyBytes Size<sup>*</sup> | ShortSeq Size<sup>*</sup> | % Reduction |
|-----------------|----------------------------|--------------------------|---------------------------|-------------|
| 0-32 nt         | 56-88 bytes                | 40-72 bytes              | 32 bytes (fixed)          | **20-64%**  |
| 33-64 nt        | 88-120 bytes               | 72-104 bytes             | 48 bytes (fixed)          | **33-60%**  |
| 65-1024 nt      | 120-1080 bytes             | 104-1064 bytes           | 48-288 bytes              | **53-73%**  |

<sup>*</sup> Object sizes were measured on Python 3.10 using `asizeof()` from the `pympler` package.

Using Cython we can move this memory representation out of Python space and into C space for better efficiency while still fulfilling the duties of a Python object. A certain amount of space goes to the Python object header and padding to 16 byte multiples for proper memory alignment, so the actual memory savings are less than the theoretical 75% reduction.

### Overhead
The time it takes to convert a list of sequences as Python bytes objects (each bytes object representing the ASCII encoding of the complete sequence) to a list of ShortSeq objects is roughly the same as it would take to do the same but with an intermediate conversion to PyUnicode. The time it takes to count the number of occurrences of each sequence in each of those lists is also roughly the same. Usually there is 5-15% runtime overhead but the current implementation is sensitive to CPU cache residency, so its range is more variable than that of PyUnicode.

### Decoding
The original ASCII form of each ShortSeq is decoded "lazily", i.e. it happens only when asked and the result isn't cached in the object, so it must be recomputed with each request. However, ShortSeqs still retain the original string's length and their ability to be compared for equality. The packed integer form of each sequence additionally acts as its pre-computed hash value. These two properties make ShortSeqs ideal for use in datastructures like sets and dictionaries, making them useful for efficiently gathering statistics involving very large volumes of short sequences.

### Longer Sequences, Featuring SIMD Encoding
Longer sequences (65 - 1024 bases) are now experimentally supported. The encoding phase uses SIMD instructions to convert 8 bases at a time for much higher throughput. Most Intel CPUs 2014 and newer, and AMD CPUs 2020 and newer support these instructions (BMI2). As a result, the PyBytes → ShortSeq conversion is often faster than PyBytes → PyUnicode. However, sequences of this length are not checked for invalid base characters, so it is the user's responsibility to ensure that only valid DNA strings are used. 
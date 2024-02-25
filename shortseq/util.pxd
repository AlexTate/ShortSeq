from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libc.stddef cimport size_t
from libc.string cimport strlen
from libc.math cimport ceil
from libcpp.cast cimport reinterpret_cast

from cpython.mem cimport PyObject_Calloc, PyObject_Free
from cpython.object cimport Py_SIZE, PyObject
from cpython.ref cimport Py_XDECREF, Py_XINCREF
from cpython.slice cimport PySlice_GetIndicesEx, PySlice_AdjustIndices
from cpython.unicode cimport PyUnicode_DecodeASCII

# For Cython, this is necessary when using these types in brackets (reinterpret_cast)
ctypedef uint64_t* llstr
ctypedef uint32_t* istr
ctypedef char* cstr

# Mask values for pext instructions
cdef uint64_t pext_mask_64
cdef uint32_t pext_mask_32

# Constants
cdef size_t NT_PER_BLOCK
cdef char[4] charmap

"""
This is a bloom filter that allows us to efficiently check
that each character is a standard base. This approach is
FAR more efficient than using set operations.
"""
cdef uint64_t bloom

"""
This is used for encoding ASCII nucleotide values to 2-bit form.
Indices correspond to the ASCII offset, and table values represent the
encoded value. While table lookup seems like it would be slower than pure 
bitwise ops, this consistently performed better (probably because the useful
portion of the table ends up in a cache line). Supports A,T,U,G,C.
"""
cdef uint8_t[91] table_91

"""
Useful for quick access to ob_sval without interacting with the CPython API
"""
ctypedef struct PyBytesObject:
    PyObject ob_base
    Py_ssize_t ob_size
    Py_hash_t ob_shash
    char ob_sval[1]

"""
For SIMD operations. 
"""
cdef extern from "x86intrin.h" nogil:
    # SSE4.2
    uint64_t _popcnt64(uint64_t __X)

    # BMI2
    uint64_t _pext_u64 (uint64_t __X, uint64_t __Y)
    uint32_t _pext_u32 (uint32_t __X, uint32_t __Y)
    uint64_t _bzhi_u64(uint64_t __X, uint32_t __Y)

"""
A little bit of hackery to allow fast access to the packed hash field of both
ShortSeq64 and ShortSeq128, since inheritance and virtual function emulation
come with a heavy cost in Cython.
"""
ctypedef struct ShortSeqGeneric:
    PyObject ob_base
    uint64_t _packed


cpdef inline printbin(header, value, value_bitwidth, chunk_bitwidth):
    """Convenience function for printing any size integer as zero-padded
    binary tokenized into chunk_bitwidth length chunks.
    
    Args:
        header: Prefix string (helps distinguish between values)
        value: The integer value to print
        value_bitwidth: The width of the value type (zero-padded to this length)
        chunk_bitwidth: The desired width of the tokenized chunks
    """
    string = f"{value:0{value_bitwidth}b}"
    chunks = [string[i:i + chunk_bitwidth] for i in range(0, len(string), chunk_bitwidth)]
    print(header + " ".join(chunks))


"""
These functions use a bloom filter to quickly check if a given
character in the provided sequence is a standard base.
NOTE: only uppercase bases supported, lowercase bases
are treated as non-standard and will be rejected.

This is a temporary fix until I find a better solution
(or perhaps this will be removed and left to user's responsibility).
"""

cdef inline bint is_base(uint8_t char) nogil:
    return bloom & (1ULL << (char & 63)) == 0


"""Validates ASCII bases 4 at a time"""

cdef inline bint _bloom_filter_32(uint32_t block) noexcept nogil:
    cdef uint32_t shifts = block & 0x3F3F3F3FL
    cdef uint64_t query = ((1ULL << ((shifts >> 0)  & 0xFF)) |
                           (1ULL << ((shifts >> 8)  & 0xFF)) |
                           (1ULL << ((shifts >> 16) & 0xFF)) |
                           (1ULL << ((shifts >> 24) & 0xFF)))

    return (bloom & query) == 0


"""Validates ASCII bases 8 at a time"""

cdef inline uint64_t _bloom_filter_64(uint64_t block) noexcept nogil:  # nolint
    cdef uint64_t shifts = block & 0x3F3F3F3F3F3F3F3FLL
    cdef uint64_t query = ((1ULL << ((shifts >> 0)  & 0xFF)) |
                           (1ULL << ((shifts >> 8)  & 0xFF)) |
                           (1ULL << ((shifts >> 16) & 0xFF)) |
                           (1ULL << ((shifts >> 24) & 0xFF)) |
                           (1ULL << ((shifts >> 32) & 0xFF)) |
                           (1ULL << ((shifts >> 40) & 0xFF)) |
                           (1ULL << ((shifts >> 48) & 0xFF)) |
                           (1ULL << ((shifts >> 56) & 0xFF)))

    return (bloom & query) == 0


"""Performs element-wise equality check for two dynamic C arrays of uint64_t's"""
cdef inline bint is_array_equal(uint64_t* a, uint64_t* b, size_t length) nogil:
    cdef size_t i

    for i in range(length):
        if a[i] != b[i]: return False

    return True


"""
Performs euclidean division and remainder and returns the result as a pair.
This is essentially the C version of Python's divmod function. Note that
Python-style modulo is performed and not C-style remainder, so this should
still give expected values for negative numbers.
"""

cdef (size_t, size_t) _divmod(size_t dividend, size_t divisor)


"""Returns the number of 64-bit blocks needed to store the specified length of bits."""
cdef size_t _bit_len_to_block_num(size_t length) noexcept nogil


"""Returns the number of 64-bit blocks needed to store the specified length of nucleotides."""
cdef size_t _nt_len_to_block_num(size_t length) noexcept nogil

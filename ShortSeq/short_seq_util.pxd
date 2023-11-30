from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libc.stddef cimport size_t
from libc.string cimport strlen
from libcpp.cast cimport reinterpret_cast

from cpython.object cimport Py_SIZE, PyObject
from cpython.ref cimport Py_XDECREF, Py_XINCREF
from cpython.unicode cimport PyUnicode_DecodeASCII


# Constants
cdef uint8_t mask
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
For vectorized operations. 
"""
cdef extern from "x86intrin.h" nogil:
    uint64_t _pext_u64 (uint64_t __X, uint64_t __Y)
    uint32_t _pext_u32 (uint32_t __X, uint32_t __Y)

"""
A little bit of hackery to allow fast access to the packed hash field of both
ShortSeq64 and ShortSeq128, since inheritance and virtual function emulation
come with a heavy cost in Cython.
"""
ctypedef struct ShortSeqGeneric:  # Todo: rename (to ShortSeqGeneric?)
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
This function uses a bloom filter to quickly check if a given
character in the provided sequence is a standard base.
NOTE: only uppercase bases supported, lowercase bases
are treated as non-standard and will be rejected.

This is a temporary fix until I find a better solution
(or perhaps this will be removed and left to user's responsibility).
"""

cdef inline bint is_base(uint8_t char) nogil:
    return bloom & (1 << (char & 63)) == 0

cdef inline bint is_array_equal(uint64_t* a, uint64_t* b, size_t length) nogil:
    cdef size_t i

    for i in range(length):
        if a[i] != b[i]: return False

    return True

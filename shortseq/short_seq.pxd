from cpython.bytes cimport PyBytes_AsString, PyBytes_Check
from cpython.unicode cimport PyUnicode_DATA, PyUnicode_GET_LENGTH, PyUnicode_Check

from libc.string cimport memcpy

from .short_seq_var cimport *
from .short_seq_128 cimport *
from .short_seq_64 cimport *
from .util cimport *

"""
Singleton values
"""

cdef ShortSeq64 empty


"""
Importable constructor functions for Cython space
"""

cpdef pack(object seq)
cdef object _from_py_str(str seq_str)
cdef object _from_py_bytes(bytes seq_bytes)
cdef object _from_chars(char* sequence)
cdef object _new(char* sequence, size_t length)

cdef ShortSeq64 _subscript(uint64_t packed, size_t offset)
cdef object _slice(uint64_t* packed, size_t offset, size_t slice_len)
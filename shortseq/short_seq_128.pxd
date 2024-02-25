from .util cimport *
from .short_seq cimport *

# Constants
cdef size_t MIN_128_NT
cdef size_t MAX_128_NT

# Cython needs help recognizing 128bit integers
cdef extern from *:
    ctypedef unsigned long long uint128_t "__uint128_t"

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_64[64]

cdef class ShortSeq128:                # 16 bytes (PyObject_HEAD)
    cdef uint128_t _packed             # 16 bytes
    cdef uint8_t _length               # 1 byte
                                       # Total: 48 bytes (15 bytes pad/free)

cdef uint128_t _marshall_bytes_128(uint8_t* seq_bytes, uint8_t length) nogil

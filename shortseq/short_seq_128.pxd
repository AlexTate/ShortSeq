from .util cimport *
from .short_seq cimport *

# Constants
cdef size_t MIN_128_NT
cdef size_t MAX_128_NT

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_64[64]

# Todo: ShortSeq128 -> ShortSeq192: 32 more bases for free!
cdef class ShortSeq128:                # 16 bytes (PyObject_HEAD)
    cdef uint64_t _packed[2]           # 16 bytes
    cdef uint8_t _length               # 1 byte
                                       # Total: 40 bytes (7 bytes pad/free)

cdef void _marshall_bytes_128(ShortSeq128 out, uint8_t* seq_bytes, uint8_t length) nogil

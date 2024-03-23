from .util cimport *
from .short_seq cimport *

# Constants
cdef size_t MIN_192_NT
cdef size_t MAX_192_NT

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_96[96]

cdef class ShortSeq192:                # 16 bytes (PyObject_HEAD)
    cdef uint64_t _packed[3]           # 24 bytes
    cdef uint8_t _length               # 1 byte
                                       # Total: 48 bytes (7 bytes pad/free)

cdef void _marshall_bytes_192(ShortSeq192 out, uint8_t* seq_bytes, uint8_t length) nogil

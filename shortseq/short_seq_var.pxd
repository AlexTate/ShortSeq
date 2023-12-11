from .util cimport *
from cpython.mem cimport PyObject_Calloc, PyObject_Free
from libc.math cimport ceil
from libc.string cimport memcmp

# Constants
cdef size_t MAX_VAR_NT
cdef size_t MIN_VAR_NT
cdef size_t MAX_REPR_LEN

# For Cython, this is necessary for using these types in brackets (reinterpret_cast)
ctypedef uint64_t* llstr
ctypedef uint32_t* istr

# Reusable buffer for unmarshalling
cdef char out_ascii_buffer_var[1024]

cdef class ShortSeqVar:                # 16 bytes (PyObject_HEAD)
    cdef uint64_t* _packed             # 8 bytes ptr
    cdef size_t _length                # 4 bytes
                                       # Total: 32 bytes
                                       #        + Heap allocation

cdef uint64_t* _marshall_bytes_var(uint8_t* seq_bytes, size_t length)
cdef uint64_t _marshall_bytes_pext_u64(uint64_t block, uint8_t* &seq_bytes, size_t n_pext) nogil
cdef uint64_t _marshall_bytes_pext_u32(uint64_t block, uint8_t* &seq_bytes, size_t n_pext) nogil
cdef uint64_t _marshall_bytes_serial(uint64_t block, uint8_t* &seq_bytes, size_t length) nogil
cdef unicode _unmarshall_bytes_var(uint64_t* enc_seq, size_t length, size_t start_block=*, size_t offset=*)

"""Validates ASCII bases 4 at a time allowing only uppercase ATGC"""

cdef inline bint _bloom_filter_32(uint32_t block) nogil:
    cdef uint32_t shifts = block & 0x3F3F3F3FL
    cdef uint64_t query = ((1UL << ((shifts >> 0) & 0xFF))  |
                           (1UL << ((shifts >> 8) & 0xFF))  |
                           (1UL << ((shifts >> 16) & 0xFF)) |
                           (1UL << ((shifts >> 24) & 0xFF)))

    return (bloom & query) == 0L

"""Validates ASCII bases 8 at a time allowing only uppercase ATGC"""

cdef inline uint64_t _bloom_filter_64(uint64_t block) nogil:
    cdef uint64_t shifts = block & 0x3F3F3F3F3F3F3F3FLL
    cdef uint64_t query = ((1ULL << ((shifts >> 0) & 0xFF))  |
                           (1ULL << ((shifts >> 8) & 0xFF))  |
                           (1ULL << ((shifts >> 16) & 0xFF)) |
                           (1ULL << ((shifts >> 24) & 0xFF)) |
                           (1ULL << ((shifts >> 32) & 0xFF)) |
                           (1ULL << ((shifts >> 40) & 0xFF)) |
                           (1ULL << ((shifts >> 48) & 0xFF)) |
                           (1ULL << ((shifts >> 56) & 0xFF)))

    return (bloom & query) == 0LL

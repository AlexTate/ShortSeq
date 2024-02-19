from .util cimport *
from cpython.mem cimport PyObject_Calloc, PyObject_Free
from libc.math cimport ceil
from libc.string cimport memcmp

# Constants
cdef size_t MAX_VAR_NT
cdef size_t MIN_VAR_NT
cdef size_t MAX_REPR_LEN

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


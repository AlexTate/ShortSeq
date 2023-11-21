from typing import Union

"""
short_seq.pxd
"""

cdef class ShortSeq:

    # No default ctor because we want to return the smallest variant possible

    def __init__(self, *args): ... # Exception?

    @staticmethod
    def from_str(str sequence): ...

    @staticmethod
    def from_bytes(bytes sequence): ...

    @staticmethod
    cdef _from_py_str(str seq_str): ...

    @staticmethod
    cdef object _from_chars(char* sequence): ...

    @staticmethod
    cdef object _from_py_bytes(bytes sequence): ...


# cdef class ShortSeqFactory:
#     # Todo: merge into ShortSeq
#
#     @staticmethod
#     cdef inline object from_bytes(bytes seqbytes): ...
#
#     @staticmethod
#     cdef inline object from_chars(char* sequence): ...
#
#     @staticmethod
#     cdef inline object _factory(char* sequence, uint8_t length): ...


cdef class ShortSeqCounter(dict):
    # Todo: get rid of count_items_py() and expand on dispatch in __(c)init__()

    def __init__(self, source): ...

    cdef _count_short_seqs(self, vector[PyObject *] &short_seqs): ...

    cdef _count_items_cc(self, vector[char *] & raw_lines): ...

    cpdef count_items_py(self, list it): ...

"""
short_seq_64.pxd
"""

cdef class ShortSeq64:
    cdef uint64_t _packed
    cdef uint8_t _length

cdef uint64_t _marshall_bytes_64(uint8_t* sequence, uint8_t length) nogil: ...
cdef inline unicode _unmarshall_bytes_64(uint64_t enc_seq, size_t length): ...

"""
short_seq_64.pxd
"""

cdef class ShortSeq128:
    cdef uint128_t _packed
    cdef uint8_t _length

cdef uint128_t _marshall_bytes_128(uint8_t* seq_bytes, uint8_t length, bint with_length=*) nogil: ...
cdef inline unicode _unmarshall_bytes_128(uint128_t enc_seq, uint8_t length = 0, bint with_length = False): ...

"""
short_seq_var.pxd
"""

cdef class ShortSeqVar:  # Todo
    cdef uint64_t *seq
    cdef uint8_t seq_len

"""
short_seq_util.pxd
"""

ctypedef struct PyBytesObject
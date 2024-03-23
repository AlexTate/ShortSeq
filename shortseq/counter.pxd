from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.exc cimport PyErr_Occurred
from cpython.list cimport PyList_GET_ITEM
from cpython.long cimport PyLong_FromSize_t

from libcpp.vector cimport vector

from .short_seq cimport _from_chars, _from_py_bytes
from .short_seq_64 cimport ShortSeq64
from .short_seq_192 cimport ShortSeq192
from .short_seq_var cimport ShortSeqVar
from .fast_read cimport *
from .util cimport *

# Singleton Values
cdef object one

cdef class ShortSeqCounter(dict):
    cdef _count_short_seq_vector(self, vector[PyObject *])
    cdef _count_chars_vector(self, vector[char*] it)
    cdef _count_py_bytes_list(self, list it)
    cdef _count_sequence(self, object seq)


cpdef ShortSeqCounter read_and_count_fastq(object filename)

"""
Private dictionary fast-path methods not currently offered by the Cython wrapper
"""

cdef extern from "Python.h":
    dict _PyDict_NewPresized(int minused)
    PyObject* _PyDict_GetItem_KnownHash(object mp, object key, Py_hash_t hash)
    PyObject* _PyDict_Pop_KnownHash(object mp, object key, Py_hash_t hash, object deflt)
    int _PyDict_SetItem_KnownHash(object mp, object key, object item, Py_hash_t hash)
    int _PyDict_DelItem_KnownHash(object mp, object key, Py_hash_t hash)
    bint _PyDict_Contains_KnownHash(object mp, object key, Py_hash_t hash)
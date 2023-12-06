# cython: language_level = 3, language=c++, profile=False, linetrace=False
import time

import cython  # For function decorators
from cython.operator cimport dereference as deref

# Singleton Values
cdef ShortSeq64 empty = ShortSeq64.__new__(ShortSeq64)
cdef object one = PyLong_FromSize_t(1)

cdef class ShortSeq:
    """Factory class that constructs the minimum object (ShortSeq64, ShortSeq128, or ShortSeqVar) for the sequence."""

    def __init__(self, *_):
        raise NotImplementedError("ShortSeq objects must be constructed using .from_* methods.")

    # === Python constructors ===============================================================

    @staticmethod
    def from_str(str seq_str):
        if not seq_str: return empty
        else: return ShortSeq._from_py_str(seq_str)

    @staticmethod
    def from_bytes(bytes seq_bytes):
        if not seq_bytes: return empty
        else: return ShortSeq._from_py_bytes(seq_bytes)

    # === Cython constructors ===============================================================

    @staticmethod
    cdef object _from_py_str(str seq_str):
        cdef bytes seq_bytes = PyUnicode_AsASCIIString(seq_str)
        return ShortSeq._from_py_bytes(seq_bytes)

    @staticmethod
    cdef object _from_py_bytes(bytes seq_bytes):
        cdef char* sequence = PyBytes_AsString(seq_bytes)
        cdef size_t length = Py_SIZE(seq_bytes)             # Note: already called PyBytes_Check() on prev line
        return ShortSeq._new(sequence, length)

    @staticmethod
    cdef inline object _from_chars(char* sequence):
        cdef size_t length = strlen(sequence) - 1           # Note: sequence is expected to have trailing \n
        return ShortSeq._new(sequence, length)

    @staticmethod
    cdef inline object _new(char* sequence, size_t length):
        if length == 0:
            return empty
        elif length <= MAX_64_NT:
            out64 = ShortSeq64.__new__(ShortSeq64)
            length = <uint8_t> length
            (<ShortSeq64> out64)._packed = _marshall_bytes_64(<uint8_t *> sequence, length)
            (<ShortSeq64> out64)._length = length
            return out64
        elif length <= MAX_128_NT:
            out128 = ShortSeq128.__new__(ShortSeq128)
            length = <uint8_t> length
            (<ShortSeq128> out128)._packed = _marshall_bytes_128(<uint8_t *> sequence, length)
            (<ShortSeq128> out128)._length = length
            return out128
        elif length <= MAX_VAR_NT:
            outvar = ShortSeqVar.__new__(ShortSeqVar)
            (<ShortSeqVar> outvar)._packed = _marshall_bytes_var(<uint8_t *> sequence, length)
            (<ShortSeqVar> outvar)._length = length
            return outvar
        else:
            raise Exception("Sequences longer than 1024 bases are not supported.")


cdef class ShortSeqCounter(dict):
    def __init__(self, source=None):
        super().__init__()

        if type(source) is list:
            self._count_py_bytes_list(source)

    def __setitem__(self, key, val):
        if type(key) not in (ShortSeq64, ShortSeq128, ShortSeqVar):
            raise TypeError(f"{self.__class__} does not support {type(key)} keys")
        PyDict_SetItem(self, key, val)

    @cython.boundscheck(False)
    cdef _count_py_bytes_list(self, list it):
        cdef bytes seqbytes

        for i in range(len(it)):
            seqbytes = <bytes> PyList_GET_ITEM(it, i)
            seq = ShortSeq._from_py_bytes(seqbytes)
            self._count_sequence(seq)

    cdef _count_short_seq_vector(self, vector[PyObject *] &short_seqs):
        for seq in short_seqs:
            self._count_sequence(<object>seq)

    @cython.boundscheck(False)
    cdef _count_chars_vector(self, vector[char *] &raw_lines):
        for seqchars in raw_lines:
            seq = ShortSeq._from_chars(seqchars)
            self._count_sequence(seq)

    cdef inline _count_sequence(self, object seq):
        cdef PyObject *oldval

        seqhash = deref(<ShortSeqGeneric*>seq)._packed
        oldval = _PyDict_GetItem_KnownHash(self, seq, seqhash)

        if oldval == NULL:
            if PyErr_Occurred():
                raise Exception("Something went wrong while retrieving sequence count.")
            if _PyDict_SetItem_KnownHash(self, seq, one, seqhash) < 0:
                raise Exception("Something went wrong while setting a new sequence count.")
        else:
            if _PyDict_SetItem_KnownHash(self, seq, <object>oldval + 1, seqhash) < 0:
                raise Exception("Something went wrong while setting an incremented sequence count.")


cpdef ShortSeqCounter read_and_count_fastq(object filename):
    cdef ShortSeqCounter counts = ShortSeqCounter()
    cdef vector[PyObject *] seqs
    # cdef vector[char *] seqs

    t1 = time.time()
    _read_fastq_short_seqs(filename.encode('utf-8'), seqs)
    # _read_fastq_chars(filename.encode('utf-8'), seqs)
    t2 = time.time()
    counts._count_short_seq_vector(seqs)
    # counts._count_chars_vector(seqs)
    t3 = time.time()

    print(f"{t2-t1:.2f}s to read {seqs.size()} total seqs, and {t3 - t2:.2f}s to count {len(counts)} unique sequences")
    return counts
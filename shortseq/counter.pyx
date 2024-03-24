import cython
import time

from cython.operator cimport dereference as deref

# Singleton Values
one = PyLong_FromSize_t(1)


cdef class ShortSeqCounter(dict):
    def __init__(self, source=None):
        super().__init__()

        if type(source) is list:
            self._count_py_bytes_list(source)

    def __setitem__(self, key, val):
        if type(key) not in (ShortSeq64, ShortSeq192, ShortSeqVar):
            raise TypeError(f"{self.__class__} does not support {type(key)} keys")
        PyDict_SetItem(self, key, val)

    @cython.boundscheck(False)
    cdef _count_py_bytes_list(self, list it):
        cdef bytes seqbytes

        for i in range(len(it)):
            seqbytes = <bytes> PyList_GET_ITEM(it, i)
            seq = _from_py_bytes(seqbytes)
            self._count_sequence(seq)

    cdef _count_short_seq_vector(self, vector[PyObject *] &short_seqs):
        for seq in short_seqs:
            self._count_sequence(<object>seq)

    @cython.boundscheck(False)
    cdef _count_chars_vector(self, vector[char *] &raw_lines):
        for seqchars in raw_lines:
            seq = _from_chars(seqchars)
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
# cython: language_level = 3, language=c++, profile=False, linetrace=False

import cython  # For function decorators
from cython.operator cimport dereference as deref

# Singleton Values
cdef ShortSeq64 empty = ShortSeq64.__new__(ShortSeq64)
# Generate 1-, 2-, and 3-mer singletons table
# kmers = ["".join(x) for n in range(1, 4) for x in itertools.product(*(["ATGC"] * n))]

# === Python constructors ===============================================================

def pack(object seq):
    if type(seq) is str:
        if not seq: return empty
        else: return _from_py_str(seq)
    elif type(seq) is bytes:
        if not seq: return empty
        else: return _from_py_bytes(seq)
    elif type(seq) is ShortSeq64:
        return seq
    elif type(seq) is ShortSeq128:
        return seq
    elif type(seq) is ShortSeqVar:
        return seq
    else:
        raise TypeError(f"Cannot pack {type(seq)} into a ShortSeq.")

def from_str(str seq_str):
    if not seq_str: return empty
    else: return _from_py_str(seq_str)

def from_bytes(bytes seq_bytes):
    if not seq_bytes: return empty
    else: return _from_py_bytes(seq_bytes)

# === Cython constructors ===============================================================

cdef object _from_py_str(str seq_str):
    cdef bytes seq_bytes = PyUnicode_AsASCIIString(seq_str)
    return _from_py_bytes(seq_bytes)

cdef object _from_py_bytes(bytes seq_bytes):
    cdef char* sequence = PyBytes_AsString(seq_bytes)
    cdef size_t length = Py_SIZE(seq_bytes)             # Note: already called PyBytes_Check() on prev line
    return _new(sequence, length)

cdef inline object _from_chars(char* sequence):
    cdef size_t length = strlen(sequence) - 1           # Note: sequence is expected to have trailing \n
    return _new(sequence, length)

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
        raise Exception(f"Sequences longer than {MAX_VAR_NT} bases are not supported.")

cdef object _from_sliced_bytes(uint64_t* packed, size_t offset, size_t length):
    """Called only when slices of ShortSeqVar or ShortSeq128 are of sufficient
    length to require a class downgrade. Slices that result in the same object
    type (as determined by length) are constructed class-local instead of here."""

    pass

cdef inline ShortSeq64 _subscript(uint64_t* packed, size_t index):
    cdef ShortSeq64 out = ShortSeq64.__new__(ShortSeq64)
    out._packed = _bzhi_u64(deref(packed) >> (index * 2), 3)
    out._length = 1
    return out

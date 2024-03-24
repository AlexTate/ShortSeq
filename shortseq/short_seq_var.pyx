# cython: language_level = 3, language=c++, profile=False, linetrace=False

cimport cython

from cython.operator cimport dereference as deref

# Importable constants
MIN_VAR_NT = 97
MAX_VAR_NT = 1024
MAX_REPR_LEN = 75

"""Used to export these constants to Python space"""
def get_domain_var(): return MIN_VAR_NT, MAX_VAR_NT

cdef class ShortSeqVar:
    def __hash__(self):
        return deref(self._packed)

    def __len__(self):
        return self._length

    def __eq__(self, other):
        if type(other) is ShortSeqVar:
            bytes_len = _nt_len_to_block_num(self._length) * sizeof(uint64_t)
            other_len = (<ShortSeqVar>other)._length
            other_ptr = (<ShortSeqVar>other)._packed
            return self._length == other_len and \
                memcmp(self._packed, <void *>other_ptr, bytes_len) == 0
        elif isinstance(other, (str, bytes)):
            return self._length == len(other) and \
                str(self) == other
        else:
            return False

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, item):
        cdef Py_ssize_t index, start, stop, step, slice_len

        if isinstance(item, slice):
            if PySlice_GetIndicesEx(item, self._length, &start, &stop, &step, &slice_len) < 0:
                raise Exception("Slice error")
            if step != 1:
                raise TypeError("Slice step not supported")
            if slice_len == 0:
                return empty
            if slice_len == 1:
                return _subscript_var(self._packed, start)

            return _slice_var(self._packed, start, slice_len)
        elif isinstance(item, int):
            index = item
            if index < 0: index += self._length
            if index < 0 or index >= self._length:
                raise IndexError("Sequence index out of range")

            return _subscript_var(self._packed, index)
        else:
            raise TypeError(f"Invalid index type: {type(item)}")

    def __str__(self):
        return _unmarshall_bytes_var(self._packed, self._length)

    def __xor__(self, ShortSeqVar other):
        if self._length != other._length:
            raise Exception(f"Hamming distance requires sequences of equal length "
                            f"({self._length} != {other._length})")

        cdef uint64_t block, block_other, block_comp
        cdef size_t n_blocks = _nt_len_to_block_num(self._length)
        cdef size_t pop_cnt = 0
        cdef size_t i

        for i in range(n_blocks):
            block = self._packed[i]
            block_other = other._packed[i]
            block_comp = block ^ block_other
            block_comp = ((block_comp >> 1) | block_comp) & 0x5555555555555555LL  # Some bases XOR to 0x3; collapse these inplace to 0x1
            pop_cnt += _popcnt64(block_comp)

        return pop_cnt

    def __sizeof__(self):
        return sizeof(ShortSeqVar) + _nt_len_to_block_num(self._length) * sizeof(uint64_t)

    def __repr__(self):
        # Truncates the sequence to MAX_REPR_LEN characters to avoid overwhelming the debugger
        cdef unicode trunc_seq = _unmarshall_bytes_var(self._packed, MAX_REPR_LEN)
        return f"<ShortSeqVar ({self._length} nt): {trunc_seq} ... >"

    def __dealloc__(self):
        if self._packed is not NULL:
            PyObject_Free(<void *>self._packed)


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_var(uint64_t* enc_seq, size_t length, size_t start_block=0, size_t offset=0):
    cdef:
        uint64_t block
        size_t blocks_spanned = _nt_len_to_block_num(length) + 1
        size_t hi = min(length, (64 - offset) // 2)
        size_t rem = length
        size_t lo = 0
        size_t i, j

    for i in range(start_block, start_block + blocks_spanned):
        block = enc_seq[i]
        if i == start_block and offset:
            block >>= offset

        for j in range(lo, hi):
            out_ascii_buffer_var[j] = charmap[block & 0b11]
            block >>= 2

        rem -= hi - lo
        lo = hi
        hi += min(rem, NT_PER_BLOCK)

    return PyUnicode_DecodeASCII(out_ascii_buffer_var, length, NULL)


cdef uint64_t* _marshall_bytes_var(uint8_t* seq_bytes, size_t length):
    cdef:
        size_t n_blocks = _nt_len_to_block_num(length)
        uint64_t* hash_arr = <llstr>PyObject_Calloc(n_blocks, sizeof(uint64_t))

    if hash_arr is NULL:
        raise MemoryError(f"Error while allocating new ShortSeq of length {length}.")

    _marshall_bytes_array(hash_arr, seq_bytes, length)
    return hash_arr


cdef inline ShortSeq64 _subscript_var(uint64_t* enc_seq, size_t index):
    """Returns a ShortSeq64 representing the specified base from the encoded sequence."""

    cdef size_t block_idx, block_offset
    block_idx, block_offset = _locate_idx(index)
    return _subscript(enc_seq[block_idx], block_offset)


cdef inline object _slice_var(uint64_t* enc_seq, size_t start, size_t slice_len):
    """Returns a new ShortSeq object representing a slice of the encoded sequence."""

    cdef size_t block_idx, block_offset
    block_idx, block_offset = _locate_idx(start)
    return _slice(enc_seq + block_idx, block_offset, slice_len)

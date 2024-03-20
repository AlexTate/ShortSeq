# cython: language_level=3, language=c++, profile=False, linetrace=False

cimport cython

""" MEASURED ON PYTHON 3.10
ShortSeq128: packs sequences up to 64 bases in length
into fixed-size objects using 2-bit encoding. For sequences
32 bases or less, ShortSeq64 is more efficient.

Max compression (64 base sequence):
    Memory footprint: 40 bytes
    PyBytes equivalent: 104 bytes (62% reduction)
    PyUnicode equivalent: 120 bytes (67% reduction)
Min compression (33 base sequence):
    Memory footprint: 40 bytes
    PyBytes equivalent: 72 bytes (44% reduction)
    PyUnicode Equivalent: 88 bytes (55% reduction)

Alternatively,
Consider encoding length into lower 6 bits (representing up to 63):
    Max sequence length: 61 bases
    Memory footprint: 40 bytes
    PyBytes equivalent: 96 bytes (58% reduction)
    PyUnicode equivalent: 112 bytes (64% reduction)
"""

# Constants
MIN_128_NT = 33
MAX_128_NT = 64

"""Used to export these constants to Python space"""
def get_domain_128(): return MIN_128_NT, MAX_128_NT

cdef class ShortSeq128:

    def __hash__(self):
        return self._packed[0]

    def __len__(self):
        return self._length

    def __eq__(self, other):
        if type(other) is ShortSeq128:
            other_len = (<ShortSeq128>other)._length
            other_ptr = (<ShortSeq128>other)._packed
            return self._length == other_len and \
                memcmp(self._packed, <void *> other_ptr, self._length) == 0
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
            elif slice_len == 1:
                return _subscript_128(self._packed, start)

            return _slice_128(self._packed, start, slice_len)
        elif isinstance(item, int):
            index = item
            if index < 0: index += self._length
            if index < 0 or index >= self._length:
                raise IndexError("Sequence index out of range")

            return _subscript_128(self._packed, index)
        else:
            raise TypeError(f"Invalid index type: {type(item)}")

    def __xor__(self, ShortSeq128 other):
        if self._length != other._length:
            raise Exception(f"Hamming distance requires sequences of equal length "
                            f"({self._length} != {other._length})")

        cdef uint64_t block, block_other, block_comp
        cdef size_t pop_cnt = 0
        cdef size_t i

        for i in range(2):
            block = self._packed[i]
            block_other = other._packed[i]
            block_comp = block ^ block_other
            block_comp = ((block_comp >> 1) | block_comp) & 0x5555555555555555LL  # Some bases XOR to 0x3; collapse these inplace to 0x1
            pop_cnt += _popcnt64(block_comp)

        return pop_cnt

    def __str__(self):
        return _unmarshall_bytes_128(self._packed, self._length)

    def __repr__(self):
        return f"<ShortSeq128 ({self._length} nt): {self}>"

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline void _marshall_bytes_128(ShortSeq128 out, uint8_t* sequence, uint8_t length) nogil:
    """Note that unlike the other _marshall_bytes_* functions, this one does not return a value. 
    Instead, it adds the packed sequence and length to the provided ShortSeq128 instance."""

    cdef:
        uint64_t* wide_iter = reinterpret_cast[llstr](sequence)
        char* nonbase_ptr
        uint8_t i

    # First 32 bases are guaranteed
    for i in reversed(range(4)):
        seq_chunk = wide_iter[i]
        if not _bloom_filter_64(seq_chunk):
            nonbase_ptr = reinterpret_cast[cstr](&seq_chunk)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 1, NULL)}")
        out._packed[0] = (out._packed[0] << 16) | _pext_u64(seq_chunk, pext_mask_64)

    for i in reversed(range(32, length)):
        seq_char = sequence[i]
        if not is_base(seq_char):
            nonbase_ptr = reinterpret_cast[cstr](&seq_char)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 1, NULL)}")
        out._packed[1] = (out._packed[1] << 2) | table_91[seq_char]

    out._length = length


@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_128(uint64_t* enc_seq, uint8_t length):
    cdef uint64_t block
    cdef uint8_t i

    # First block is guaranteed full
    block = enc_seq[0]
    for i in range(32):
        out_ascii_buffer_64[i] = charmap[block & 0b11]
        block >>= 2

    block = enc_seq[1]
    for i in range(32, length):
        out_ascii_buffer_64[i] = charmap[block & 0b11]
        block >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_64, length, NULL)


cdef inline ShortSeq64 _subscript_128(uint64_t* enc_seq, size_t index):
    """Returns a ShortSeq64 representing the specified base from the encoded sequence."""

    cdef size_t block_idx, block_offset
    block_idx, nt_offset = _divmod(index, 32)
    return _subscript(enc_seq[block_idx], nt_offset)


cdef inline object _slice_128(uint64_t* enc_seq, size_t start, size_t slice_len):
    """Returns a new ShortSeq object representing a slice of the encoded sequence."""

    cdef size_t block_idx, offset
    block_idx, offset = _divmod(start * 2, 64)
    return _slice(enc_seq + block_idx, offset, slice_len)
# cython: language_level=3, language=c++, profile=False, linetrace=False

cimport cython

""" MEASURED ON PYTHON 3.10
ShortSeq128: packs sequences up to 64 bases in length
into fixed-size objects using 2-bit encoding. For sequences
32 bases or less, ShortSeq64 is more efficient.

Max compression (64 base sequence):
    Memory footprint: 48 bytes
    PyBytes equivalent: 104 bytes (54% reduction)
    PyUnicode equivalent: 120 bytes (60% reduction)
Min compression (33 base sequence):
    Memory footprint: 48 bytes
    PyBytes equivalent: 72 bytes (33% reduction)
    PyUnicode Equivalent: 88 bytes (45% reduction)

Alternatively,
Consider encoding length into lower 6 bits (representing up to 63):
    Max sequence length: 61 bases
    Memory footprint: 48 bytes
    PyBytes equivalent: 96 bytes (50% reduction)
    PyUnicode equivalent: 112 bytes (57% reduction)
"""

# Constants
MIN_128_NT = 33
MAX_128_NT = 64

"""Used to export these constants to Python space"""
def get_domain_128(): return MIN_128_NT, MAX_128_NT

cdef class ShortSeq128:

    def __hash__(self):
        return <uint64_t>self._packed

    def __len__(self):
        return self._length

    def __eq__(self, other):
        if type(other) is ShortSeq128:
            return self._length == (<ShortSeq128> other)._length and \
                   self._packed == (<ShortSeq128> other)._packed
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

        cdef uint128_t xor = self._packed ^ (<ShortSeq128> other)._packed
        cdef uint64_t lo = <uint64_t> xor
        cdef uint64_t hi = xor >> 64

        # Some bases XOR to 0x3; collapse these results to 0x1 inplace
        lo = ((lo >> 1) | lo) & 0x5555555555555555LL
        hi = ((hi >> 1) | hi) & 0x5555555555555555LL

        return _popcnt64(lo) + _popcnt64(hi)

    def __str__(self):
        return _unmarshall_bytes_128(self._packed, self._length)

    def __repr__(self):
        return f"<ShortSeq128 ({self._length} nt): {self}>"


@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint128_t _marshall_bytes_128(uint8_t* sequence, uint8_t length) nogil:
    """Encodes the input sequence into a 128 bit integer using 2-bit encoding.
    Note that from a high-level memory perspective, the layout of the resulting marshalled
    bytes is more similar to that of ShortSeq64, which stores bases in a sequentially continuous
    manner, than it is to ShortSeqVar, which stores bases in uint64_t array elements."""

    cdef:
        uint64_t* wide_iter = reinterpret_cast[llstr](sequence)
        uint64_t hi, lo
        char* nonbase_ptr
        uint8_t i

    # Last 32 bases can be variable length. Handle them one at a time.
    for i in reversed(range(32, length)):
        seq_char = sequence[i]
        if not is_base(seq_char):
            nonbase_ptr = reinterpret_cast[cstr](&seq_char)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 1, NULL)}")
        hi = (hi << 2) | table_91[seq_char]

    # First 32 bases are always guaranteed. Handle them 8 at a time.
    for i in reversed(range(4)):
        seq_chunk = wide_iter[i]
        if not _bloom_filter_64(seq_chunk):
            nonbase_ptr = reinterpret_cast[cstr](&seq_chunk)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 1, NULL)}")
        lo = (lo << 16) | _pext_u64(seq_chunk, pext_mask_64)

    return (<uint128_t>hi << 64) | <uint128_t>lo

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_128(uint128_t enc_seq, uint8_t length):
    cdef uint8_t i

    for i in range(length):
        out_ascii_buffer_64[i] = charmap[enc_seq & 0b11]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_64, length, NULL)


cdef inline ShortSeq64 _subscript_128(uint128_t enc_seq, size_t index):
    """Returns a ShortSeq64 representing the specified base from the encoded sequence."""

    cdef uint64_t half
    cdef size_t corr_idx

    if index < 32:
        half = <uint64_t> enc_seq
        corr_idx = index
    else:
        half = enc_seq >> 64
        corr_idx = index - 32

    return _subscript(half, corr_idx)


cdef inline object _slice_128(uint128_t enc_seq, size_t start, size_t slice_len):
    """Returns a new ShortSeq object representing a slice of the encoded sequence."""

    cdef uint64_t* block_ptr = reinterpret_cast[llstr](&enc_seq)
    cdef size_t block_idx, offset

    block_idx, offset = _divmod(start * 2, 64)
    return _slice(block_ptr + block_idx, offset, slice_len)

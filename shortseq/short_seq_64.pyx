# cython: language_level = 3, language=c++, profile=False, linetrace=False

cimport cython

""" MEASURED ON PYTHON 3.10
ShortSeq64: packs sequences up to 32 bases in length
into fixed-size objects using 2-bit encoding.

Max compression (32 base sequence):
    Memory footprint: 32 bytes
    PyBytes equivalent: 72 bytes (56% reduction)
    PyUnicode equivalent: 88 bytes (64% reduction)
Min compression (empty sequence):
    Memory footprint: 32 bytes
    PyBytes equivalent: 40 bytes (20% reduction)
    PyUnicode equivalent: 56 bytes (42% reduction)
    
Alternatively,
Consider if length was encoded into the lower 6 bits (representing up to 63):
    Max sequence length: 29 bases
    Memory footprint: 32 bytes
    PyBytes equivalent: 64 bytes (50% reduction)
    PyUnicode equivalent: 80 bytes (60% reduction)
"""

# Constants
MIN_64_NT = 0
MAX_64_NT = 32

"""Used to export these constants to Python space"""
def get_domain_64(): return MIN_64_NT, MAX_64_NT

cdef class ShortSeq64:
    # Todo: decide whether to standardize or remove (esp. setter...)

    @property
    def _length(self) -> uint8_t:
        return self._length

    @_length.setter
    def _length(self, uint8_t val):
        self._length = val

    @property
    def _packed(self) -> uint64_t:
        return self._packed

    @_packed.setter
    def _packed(self, val):
        self._packed = val

    def __hash__(self) -> uint64_t:
        return self._packed

    def __len__(self):
        return self._length

    def __eq__(self, other):
        if type(other) is ShortSeq64:
            return self._length == (<ShortSeq64>other)._length and \
                   self._packed == (<ShortSeq64>other)._packed
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
            return _unmarshall_bytes_64(self._packed >> (start * 2), slice_len)
        elif isinstance(item, int):
            index = item
            if index < 0: index += self._length
            if index < 0 or index >= self._length:
                raise IndexError("Sequence index out of range")
            return _unmarshall_bytes_64(self._packed >> (index * 2), 1)
        else:
            raise TypeError(f"Invalid index type: {type(item)}")

    def __xor__(self, ShortSeq64 other):
        if self._length != other._length:
            raise Exception("Hamming distance requires sequences of equal length")

        cdef uint64_t comp = self._packed ^ (<ShortSeq64>other)._packed
        comp = ((comp >> 1) | comp) & 0x5555555555555555LL  # Some bases XOR to 0x3; collapse these inplace to 0x1
        return _popcnt64(comp)

    def __str__(self):
        return _unmarshall_bytes_64(self._packed, self._length)

    def __repr__(self):
        return f"<ShortSeq64 ({self._length} nt): {self}>"


@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint64_t _marshall_bytes_64(uint8_t* sequence, uint8_t length) nogil:
    cdef uint64_t hashed = 0LL
    cdef char* nonbase_ptr
    cdef uint8_t i

    for i in reversed(range(length)):
        seq_char = sequence[i]
        if not is_base(seq_char):
            nonbase_ptr = reinterpret_cast[cstr](&seq_char)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 1, NULL)}")
        hashed = (hashed << 2) | table_91[seq_char]

    return hashed

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_64(uint64_t enc_seq, size_t length):
    cdef uint8_t i

    for i in range(length):
        out_ascii_buffer_32[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_32, length, NULL)
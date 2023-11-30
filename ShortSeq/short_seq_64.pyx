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

    def __eq__(self, other):
        if type(other) is ShortSeq64:
            return self._length == (<ShortSeq64>other)._length and \
                   self._packed == (<ShortSeq64>other)._packed
        else:
            return False

    def __str__(self):
        return _unmarshall_bytes_64(self._packed, self._length)

    def __len__(self):
        return self._length

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint64_t _marshall_bytes_64(uint8_t* sequence, uint8_t length) nogil:
    cdef uint64_t hashed = 0L
    cdef uint8_t seq_char
    cdef uint8_t i

    for i in reversed(range(length)):
        seq_char = sequence[i]
        if not is_base(seq_char):
            raise Exception(f"Unsupported base character: {seq_char}")
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
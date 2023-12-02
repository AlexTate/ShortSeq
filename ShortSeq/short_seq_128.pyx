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

    def __str__(self):
        return _unmarshall_bytes_128(self._packed, self._length)

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint128_t _marshall_bytes_128(uint8_t* sequence, uint8_t length, bint with_length = False) nogil:
    cdef uint128_t hashed = 0LL
    cdef uint8_t seq_char
    cdef uint8_t i

    for i in reversed(range(length)):
        seq_char = sequence[i]
        if not is_base(seq_char):
            raise Exception(f"Unsupported base character: {seq_char}")
        hashed = (hashed << 2) | table_91[seq_char]

    if with_length:
        hashed = (hashed << 8) | length

    return hashed

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_128(uint128_t enc_seq, uint8_t length = 0, bint with_length = False):
    cdef uint8_t i

    if with_length:
        length = enc_seq & 0xFF
        enc_seq >>= 8

    for i in range(length):
        out_ascii_buffer_64[i] = charmap[enc_seq & mask]
        enc_seq >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_64, length, NULL)
# cython: language_level = 3, language=c++, profile=False, linetrace=False

cimport cython

from cython.operator cimport dereference as deref

cdef class ShortSeqVar:
    def __hash__(self):
        return deref(self._packed)

    def __len__(self):
        return self._length

    def __eq__(self, other):
        if type(other) is ShortSeqVar:
            other_len = (<ShortSeqVar>other)._length
            other_ptr = (<ShortSeqVar>other)._packed
            return self._length == other_len and \
                memcmp(self._packed, <void *>other_ptr, self._length) == 0
        elif isinstance(other, (str, bytes)):
            return self._length == len(other) and \
                   str(self) == other
        else:
            return False

    def __str__(self):
        return _unmarshall_bytes_var(self._packed, self._length)

    def __repr__(self):
        # Clips the sequence to MAX_REPR_LEN characters to avoid overwhelming the debugger
        cdef unicode clipped_seq = _unmarshall_bytes_var(self._packed, MAX_REPR_LEN)
        return f"<ShortSeqVar ({self._length} nt): {self} ... >"

    def __dealloc__(self):
        if self._packed is not NULL:
            PyObject_Free(<void *>self._packed)


cdef inline unicode _unmarshall_bytes_var(uint64_t* enc_seq, size_t length):
    cdef:
        uint64_t block
        size_t num_blocks = <size_t>ceil(<double>length / 32.0)
        size_t i, j, lo, hi

    for i in range(num_blocks):
        block = enc_seq[i]
        lo = i * 32
        hi = lo + min(32, length)
        for j in range(lo, hi):
            out_ascii_buffer_var[j] = charmap[block & mask]
            block >>= 2

    return PyUnicode_DecodeASCII(out_ascii_buffer_var, length, NULL)

@cython.wraparound(False)
cdef uint64_t* _marshall_bytes_var(uint8_t* seq_bytes, size_t length):
    cdef:
        size_t n_blocks = <size_t>ceil(<double>length / 32.0)
        uint64_t* hash_arr = <llstr>PyObject_Calloc(n_blocks, sizeof(uint64_t))
        uint8_t* seq_it = seq_bytes
        size_t seq_rem = length
        size_t block_rem, i
        uint64_t block

        size_t offset_p64, offset_p32
        size_t n_p64, p64_rem, n_p32, serial

    if hash_arr is NULL:
        raise MemoryError(f"Error while allocating new ShortSeq of length {length}.")

    for i in range(n_blocks):
        block_rem = min(32, seq_rem)
        seq_rem -= block_rem
        block = 0LL

        n_p64 = <size_t>(block_rem / 8)
        p64_rem = block_rem % 8

        n_p32 = <size_t>(p64_rem / 4)
        serial = p64_rem % 4

        offset_p64 = n_p64 * 8
        offset_p32 = n_p32 * 4

        if serial:
            ser_it = <uint8_t*>(seq_it + offset_p32 + offset_p64)
            block = _marshall_bytes_serial(block, ser_it, serial)    # should occupy upper bits, but be last bases
        if n_p32:
            p32_it = <uint8_t*>(seq_it + offset_p64)
            block = _marshall_bytes_pext_u32(block, p32_it, n_p32)
        if n_p64:
            block = _marshall_bytes_pext_u64(block, seq_it, n_p64)   # should occupy lower bits, but be first bases

        hash_arr[i] = block
        seq_it += block_rem

    return hash_arr

cdef uint64_t pext_mask_64 = 0x0606060606060606

@cython.wraparound(False)
cdef inline uint64_t _marshall_bytes_pext_u64(uint64_t block, uint8_t* &seq_bytes, size_t n_pext) nogil:
    """Packs 8 bases at a time into the uint64_t `block`, producing a maximum of 64 bits (for 32 bases)
    and a minimum of 16 bits (for 8 bases). 
    """

    cdef:
        uint64_t* sequence = reinterpret_cast[llstr](seq_bytes)
        size_t i

    for i in reversed(range(n_pext)):
        block = (block << 16) | _pext_u64(sequence[i], pext_mask_64)

    return block


cdef uint32_t pext_mask_32 = 0x06060606

@cython.wraparound(False)
cdef inline uint64_t _marshall_bytes_pext_u32(uint64_t block, uint8_t* &seq_bytes, size_t n_pext) nogil:
    """Packs 4 bases at a time into the uint64_t `block`"""

    cdef:
        uint32_t* sequence = reinterpret_cast[istr](seq_bytes)
        size_t i

    for i in reversed(range(n_pext)):
        block = (block << 8) | _pext_u32(sequence[i], pext_mask_32)

    return block

@cython.wraparound(False)
cdef inline uint64_t _marshall_bytes_serial(uint64_t block, uint8_t* &seq_bytes, size_t length) nogil:
    """Packs bases one at a time (sequential, non-vector) into the uint64_t `block`"""

    cdef uint8_t seq_char
    cdef uint8_t i

    for i in reversed(range(length)):
        seq_char = seq_bytes[i]
        if not is_base(seq_char):
            raise Exception(f"Unsupported base character: {seq_char}")
        block = (block << 2) | table_91[seq_char]

    return block
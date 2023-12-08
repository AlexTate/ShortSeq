# cython: language_level = 3, language=c++, profile=False, linetrace=False

cimport cython

from cython.operator cimport dereference as deref

# Importable constants
MIN_VAR_NT = 65
MAX_VAR_NT = 1024
MAX_REPR_LEN = 75

# Constants
cdef size_t NT_PER_BLOCK = 32

"""Used to export these constants to Python space"""
def get_domain_var(): return MIN_VAR_NT, MAX_VAR_NT

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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, item):
        cdef Py_ssize_t index, start, stop, step, slice_len
        cdef size_t block_idx, block_offset

        if isinstance(item, slice):
            if PySlice_GetIndicesEx(item, self._length, &start, &stop, &step, &slice_len) < 0:
                raise Exception("Slice error")
            if step != 1:
                raise TypeError("Slice step not supported")

            block_idx, block_offset = _locate_idx(start)
            return _unmarshall_bytes_var(self._packed, slice_len, block_idx, block_offset)
        elif isinstance(item, int):
            index = item
            if index < 0: index += self._length
            if index < 0 or index >= self._length:
                raise IndexError("Sequence index out of range")

            block_idx, block_offset = _locate_idx(index)
            return _unmarshall_bytes_var(self._packed, 1, block_idx, block_offset)
        else:
            raise TypeError(f"Invalid index type: {type(item)}")

    def __str__(self):
        return _unmarshall_bytes_var(self._packed, self._length)

    def __xor__(self, ShortSeqVar other):
        cdef uint64_t block, block_other, block_comp
        cdef size_t n_blocks = _length_to_block_num(self._length)
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
        return sizeof(ShortSeqVar) + _length_to_block_num(self._length) * sizeof(uint64_t)

    def __repr__(self):
        # Clips the sequence to MAX_REPR_LEN characters to avoid overwhelming the debugger
        cdef unicode clipped_seq = _unmarshall_bytes_var(self._packed, MAX_REPR_LEN)
        return f"<ShortSeqVar ({self._length} nt): {self} ... >"

    def __dealloc__(self):
        if self._packed is not NULL:
            PyObject_Free(<void *>self._packed)

@cython.cdivision(True)
cdef inline size_t _length_to_block_num(size_t length):
    """Returns the number of 64-bit blocks needed to store the given length."""

    return <size_t>ceil(<double>length / <double>NT_PER_BLOCK)

@cython.cdivision(True)
cdef inline (size_t, size_t) _locate_idx(size_t index):
    """Returns the block index and offset (in packed units) where the given index (in nt units) is located."""

    cdef size_t block_idx = index // NT_PER_BLOCK
    cdef size_t block_offset = 2 * (index % NT_PER_BLOCK)
    return block_idx, block_offset

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline unicode _unmarshall_bytes_var(uint64_t* enc_seq, size_t length, size_t start_block=0, size_t offset=0):
    cdef:
        uint64_t block
        size_t blocks_spanned = _length_to_block_num(length) + 1
        size_t hi = min(length, (64 - offset) // 2)
        size_t rem = length
        size_t lo = 0
        size_t i, j

    for i in range(start_block, start_block + blocks_spanned):
        block = enc_seq[i]
        if i == start_block and offset:
            block >>= offset

        for j in range(lo, hi):
            out_ascii_buffer_var[j] = charmap[block & mask]
            block >>= 2

        rem -= hi - lo
        lo = hi
        hi += min(rem, NT_PER_BLOCK)

    return PyUnicode_DecodeASCII(out_ascii_buffer_var, length, NULL)

@cython.wraparound(False)
cdef uint64_t* _marshall_bytes_var(uint8_t* seq_bytes, size_t length):
    cdef:
        size_t n_blocks = _length_to_block_num(length)
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
        block_rem = min(seq_rem, NT_PER_BLOCK)
        seq_rem -= block_rem
        block = 0LL

        n_p64, p64_rem = _divmod(block_rem, 8)
        n_p32, serial = _divmod(p64_rem, 4)
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
import cython

@cython.cdivision(True)
cdef inline (size_t, size_t) _divmod(size_t dividend, size_t divisor) nogil:
    cdef size_t div_res
    cdef size_t mod_res

    if divisor == 0: raise ZeroDivisionError()
    div_res = dividend // divisor
    mod_res = (((dividend % divisor) + divisor) % divisor)

    return div_res, mod_res

@cython.cdivision(True)
@cython.exceptval(check=False)
cdef inline (size_t, size_t) _locate_idx(size_t index) nogil:
    """Returns the block index and bit offset where the index (given in nt units) is located."""

    cdef size_t block_idx, block_offset
    block_idx, block_offset = _divmod(index, NT_PER_BLOCK)
    return block_idx, block_offset * 2

@cython.cdivision(True)
cdef inline size_t _bit_len_to_block_num(size_t length) noexcept nogil:
    """Returns the number of 64-bit blocks needed to store the specified length of bits."""

    return <size_t>ceil(<double>length / 64.0)

@cython.cdivision(True)
cdef inline size_t _nt_len_to_block_num(size_t length) noexcept nogil:
    """Returns the number of 64-bit blocks needed to store the specified length of nucleotides."""

    return <size_t>ceil(<double>length / <double>NT_PER_BLOCK)

"""
CONSTANTS
"""

cdef uint64_t pext_mask_64 = 0x0606060606060606
cdef uint32_t pext_mask_32 = 0x06060606

cdef size_t NT_PER_BLOCK = 32

cdef uint8_t[91] table_91 = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4
]

cdef char[4] charmap = [b'A', b'C', b'T', b'G']

""" To create the bloom filter:
lc_bases = 'atgc'
uc_bases = 'ATGC'
bases = [ord(c) for c in lc_bases + uc_bases]

lc_alpha_nb = [c for c in range(ord('a'), ord('z')+1) if chr(c) not in lc_bases]
uc_alpha_nb = [c for c in range(ord('A'), ord('Z')+1) if chr(c) not in uc_bases]

special_a = list(range(ord(' '), ord('/')+1))
special_b = list(range(ord(':'), ord('@')+1))
special_c = list(range(ord('['), ord('`')+1))
special_d = list(range(ord('{'), ord('~')+1))

alpha = lc_alpha_nb + uc_alpha_nb
numeric = list(range(ord('0'), ord('9')+1))
special = special_a + special_b + special_c + special_d

bloom = 0
for nb in alpha + numeric + special:
    bloom |= 1 << (nb & 63)
"""
cdef uint64_t bloom = 0xFFFFFFFFFFEFFF75


@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint64_t _marshall_full_block(uint8_t* sequence) nogil:
    """Encodes 32 nucleotides into a uint64_t block efficiently."""

    cdef:
        uint64_t * chunk_iter = reinterpret_cast[llstr](sequence)
        uint64_t block = 0ULL
        char * nonbase_ptr
        uint8_t i

    for i in reversed(range(4)):
        chunk = chunk_iter[i]
        if not _bloom_filter_64(chunk):
            nonbase_ptr = reinterpret_cast[cstr](&chunk)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 8, NULL)}")
        block = (block << 16) | _pext_u64(chunk, pext_mask_64)

    return block


@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline uint64_t _marshall_partial_block(uint8_t * sequence, size_t length) nogil:
    """Encodes less than 32 nucleotides into a uint64_t block."""

    cdef:
        uint64_t block = 0ULL
        char * nonbase_ptr
        uint8_t i

    for i in reversed(range(length)):
        seq_char = sequence[i]
        if not is_base(seq_char):
            nonbase_ptr = reinterpret_cast[cstr](&seq_char)
            raise Exception(f"Unsupported base character: {PyUnicode_DecodeASCII(nonbase_ptr, 1, NULL)}")
        block = (block << 2) | table_91[seq_char]

    return block

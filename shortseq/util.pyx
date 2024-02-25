import cython

@cython.cdivision(True)
cdef inline (size_t, size_t) _divmod(size_t dividend, size_t divisor):
    cdef size_t div_res
    cdef size_t mod_res

    if divisor == 0: raise ZeroDivisionError()
    div_res = dividend // divisor
    mod_res = (((dividend % divisor) + divisor) % divisor)

    return div_res, mod_res

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

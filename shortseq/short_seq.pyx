# cython: language_level = 3, language=c++, profile=False, linetrace=False

import cython  # For function decorators
from cython.operator cimport dereference as deref

# Singleton Values
cdef ShortSeq64 empty = ShortSeq64.__new__(ShortSeq64)
# Generate 1-, 2-, and 3-mer singletons table
# kmers = ["".join(x) for n in range(1, 4) for x in itertools.product(*(["ATGC"] * n))]

# === Python constructors ===============================================================

@cython.always_allow_keywords(False)
cpdef pack(object seq):
    if PyUnicode_Check(seq):
        if not seq: return empty
        else: return _from_py_str(seq)
    elif PyBytes_Check(seq):
        if not seq: return empty
        else: return _from_py_bytes(seq)
    elif type(seq) is ShortSeq64:
        return seq
    elif type(seq) is ShortSeq128:
        return seq
    elif type(seq) is ShortSeqVar:
        return seq
    else:
        raise TypeError(f'Cannot pack objects of type "{type(seq)}"')

def from_str(str seq_str):
    if not seq_str: return empty
    else: return _from_py_str(seq_str)

def from_bytes(bytes seq_bytes):
    if not seq_bytes: return empty
    else: return _from_py_bytes(seq_bytes)

# === Cython constructors ===============================================================

cdef object _from_py_str(str seq_str):
    cdef char* sequence = <char *>PyUnicode_DATA(seq_str)
    cdef size_t length = PyUnicode_GET_LENGTH(seq_str)
    return _new(sequence, length)

cdef object _from_py_bytes(bytes seq_bytes):
    cdef char* sequence = PyBytes_AsString(seq_bytes)
    cdef size_t length = Py_SIZE(seq_bytes)             # Note: already called PyBytes_Check() on prev line
    return _new(sequence, length)

cdef inline object _from_chars(char* sequence):
    cdef size_t length = strlen(sequence) - 1           # Note: sequence is expected to have trailing \n
    return _new(sequence, length)

cdef inline object _new(char* sequence, size_t length):
    if length == 0:
        return empty
    elif length <= MAX_64_NT:
        out64 = ShortSeq64.__new__(ShortSeq64)
        length = <uint8_t> length
        (<ShortSeq64> out64)._packed = _marshall_bytes_64(<uint8_t *> sequence, length)
        (<ShortSeq64> out64)._length = length
        return out64
    elif length <= MAX_128_NT:
        out128 = ShortSeq128.__new__(ShortSeq128)
        length = <uint8_t> length
        (<ShortSeq128> out128)._packed = _marshall_bytes_128(<uint8_t *> sequence, length)
        (<ShortSeq128> out128)._length = length
        return out128
    elif length <= MAX_VAR_NT:
        outvar = ShortSeqVar.__new__(ShortSeqVar)
        (<ShortSeqVar> outvar)._packed = _marshall_bytes_var(<uint8_t *> sequence, length)
        (<ShortSeqVar> outvar)._length = length
        return outvar
    else:
        raise Exception(f"Sequences longer than {MAX_VAR_NT} bases are not supported.")


cdef inline ShortSeq64 _subscript(uint64_t packed, size_t index):
    """Constructs a ShortSeq64 object from a single base of a bit-packed sequence.
    
    Args:
        packed: uint64_t
            The packed sequence to be subscripted.
        index: size_t
            The index (in nt units) of the encoded base.
    """

    cdef ShortSeq64 out = ShortSeq64.__new__(ShortSeq64)
    out._packed = (packed >> (index * 2)) & 0b11
    out._length = 1
    return out


cdef inline object _slice(uint64_t* packed, size_t offset, size_t slice_len_nts):
    """Constructs a new ShortSeq as a substring of an existing bit-packed sequence.
    
    Args:
        packed: A pointer to the block in the bit-packed sequence
            which contains the first base of the slice interval.
        offset: The offset in bits to the first base of the slice.
        slice_len_nts: The length of the slice in nucleotides.
        
    Returns:
        A new ShortSeq object containing the specified slice. The 
        type of the returned object depends on the length of the 
        slice and may not match the type of the original ShortSeq.
    """

    if slice_len_nts <= MAX_64_NT:
        return _slice_to_ShortSeq64(packed, offset, slice_len_nts)
    elif slice_len_nts <= MAX_128_NT:
        return _slice_to_ShortSeq128(packed, offset, slice_len_nts)
    elif slice_len_nts <= MAX_VAR_NT:
        return _slice_to_ShortSeqVar(packed, offset, slice_len_nts)
    else:
        raise Exception(f"Slice length {slice_len_nts} exceeds the max sequence length ({MAX_VAR_NT}).")


cdef ShortSeq64 _slice_to_ShortSeq64(uint64_t* packed, size_t offset, size_t slice_len_nts):
    """Constructs a new ShortSeq64 object from a bit-packed sequence slice.
    
    Args:
        packed: A pointer to the block in the bit-packed sequence
            which contains the first base of the slice interval.
        offset: The offset, in nucleotides, to the first base 
            of the slice (i.e. the number of bases to skip in
            the first block). Valid values are [0, 32).
        slice_len_nts: The length of the slice in nucleotides.
    """

    cdef:
        ShortSeq64 out64 = ShortSeq64.__new__(ShortSeq64)
        size_t slice_len_bits = slice_len_nts * 2
        uint64_t result

    if offset + slice_len_bits > 64:
        result = (packed[0] >> offset) | (packed[1] << (64 - offset))
        result = _bzhi_u64(result, slice_len_bits)
    else:
        result = _bzhi_u64(packed[0] >> offset, slice_len_bits)

    out64._packed = result
    out64._length = <uint8_t> slice_len_nts
    return out64


cdef ShortSeq128 _slice_to_ShortSeq128(uint64_t* packed, size_t offset, size_t slice_len_nts):
    """Constructs a new ShortSeq128 object from a bit-packed sequence slice.
    
    Args:
        packed: A pointer to the block in the bit-packed sequence
            which contains the first base of the slice interval.
        offset: The offset, in nucleotides, to the first base 
            of the slice (i.e. the number of bases to skip in
            the first block). Valid values are [0, 32).
        slice_len_nts: The length of the slice in nucleotides.
    """

    cdef:
        ShortSeq128 out128 = ShortSeq128.__new__(ShortSeq128)
        size_t slice_len_bits = slice_len_nts * 2
        cdef uint64_t* result

    result = reinterpret_cast[llstr](&out128._packed)
    _shift_copy_trim(result, packed, slice_len_bits, offset)

    out128._length = <uint8_t> slice_len_nts
    return out128

cdef ShortSeqVar _slice_to_ShortSeqVar(uint64_t* packed, size_t offset, size_t slice_len_nts):
    """Constructs a new ShortSeqVar object from a bit-packed sequence slice.
    
    Args:
        packed: A pointer to the block in the bit-packed sequence
            which contains the first base of the slice interval.
        offset: The offset, in nucleotides, to the first base 
            of the slice (i.e. the number of bases to skip in
            the first block). Valid values are [0, 32).
        slice_len_nts: The length of the slice in nucleotides.
    """

    cdef:
        ShortSeqVar outvar = ShortSeqVar.__new__(ShortSeqVar)
        size_t n_blocks = _nt_len_to_block_num(slice_len_nts)
        uint64_t* result = <llstr> PyObject_Calloc(n_blocks, sizeof(uint64_t))
        size_t slice_len_bits = slice_len_nts * 2

    if result is NULL:
        raise MemoryError(f"Error while allocating new ShortSeq of length {slice_len_nts}.")

    try:
        _shift_copy_trim(result, packed, slice_len_bits, offset)
    except Exception as e:
        PyObject_Free(result)
        raise e

    outvar._packed = result
    outvar._length = slice_len_nts
    return outvar


cdef inline void _shift_copy_trim(uint64_t* dst, uint64_t* src, size_t slice_len_bits, size_t offset=0):
    """Copies the specified number of bits from src to dst with an optional offset.
    If an offset is specified, it is applied to the entire copied length and
    blocks are reassembled to produce a valid bit-packed sequence. The final
    block of the copied sequence is trimmed to the correct length if necessary.
    
    It is important to zero out, or trim, the final block to the correct length.
    Otherwise, operations that take place on entire blocks
    
    Args:
        dst: uint64_t*
            The destination buffer to copy the sequence to.
        src: uint64_t*
            The source buffer to copy the sequence from.
        slice_len_bits: size_t
            The length of the region to copy in bits.
        offset: size_t
            The offset in bits to the first base of the slice.  
    """

    cdef:
        size_t n_blocks_dst = _bit_len_to_block_num(slice_len_bits)
        size_t tail = slice_len_bits % 64
        size_t complement = 64 - offset
        size_t i

    if offset == 0:
        # Copy directly from the source
        memcpy(dst, src, n_blocks_dst * sizeof(uint64_t))
    else:
        # Copy and reassemble blocks with offset applied to the entire length
        for i in range(n_blocks_dst):
            dst[i] = (src[i] >> offset) | (src[i + 1] << complement)

    if tail:
        # Trim the final block to the correct length
        dst[n_blocks_dst - 1] = _bzhi_u64(dst[n_blocks_dst - 1], tail)

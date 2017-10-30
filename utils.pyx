
SNP_MUTATION_CODES = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

cdef int _OTHER_MUTATION_CODE = 6
cdef int _INDEL_MUTATION_CODE = 5
cdef int _SUBST_MUTATION_CODE = 0


def make_search_key(chrom, pos, ref, alt):
    cdef int ref_length = len(ref)
    cdef int alt_length = len(alt)
    cdef long int mutation_code = _OTHER_MUTATION_CODE

    if ref_length == 1 and alt_length == 1:
        if alt in SNP_MUTATION_CODES:
            mutation_code = SNP_MUTATION_CODES[alt]
    elif ref_length != alt_length:
        mutation_code = _INDEL_MUTATION_CODE
    elif ref_length == alt_length:
        mutation_code = _SUBST_MUTATION_CODE

    # chrom - uint8 <- be careful with this part, because
    # pytables do not support uint64 as index and we cut
    # last bit when writing data (technically chrome should
    # be considered as uint7)
    # pos - uint32
    # other - (mutation_code)
    return chrom << 56 | pos << 24 | mutation_code

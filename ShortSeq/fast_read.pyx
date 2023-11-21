

cdef inline void _read_fastq_short_seqs(char* fname, vector[PyObject *] &out):
    cdef:
        FILE *cfile = fopen(fname, <char*>"rb")
        char *line = NULL
        size_t count = 1
        size_t l = 0

    if cfile == NULL:
        raise Exception(f"{str(fname)}: Something went wrong while reading this file.")

    while getline(&line, &l, cfile) != -1:
        if count % 2 == 0 and count % 4 != 0:
            seq = ShortSeq._from_chars(line)
            Py_XINCREF(<PyObject *>seq)
            out.push_back(<PyObject *>seq)

        count += 1
    fclose(cfile)


cdef inline void _read_fastq_chars(char * fname, vector[char *] &out) nogil:
    cdef:
        FILE *cfile = fopen(fname, "rb")
        char *linecpy = NULL
        char *line = NULL
        size_t count = 1
        size_t l = 0

    if cfile == NULL:
        raise Exception(f"{str(fname)}: Something went wrong while reading this file.")

    while getline(&line, &l, cfile) != -1:
        if count % 2 == 0 and count % 4 != 0:
            linecpy = strdup(line)
            out.push_back(linecpy)

        count += 1
    fclose(cfile)
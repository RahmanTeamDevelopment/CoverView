cdef extern from "sam.h":
    ctypedef struct bam1_t


cdef class ReadArray:
    cdef void append(self, bam1_t* read)
    cdef void setWindowPointers(self, int start, int end, bam1_t** window_start, bam1_t** window_end)
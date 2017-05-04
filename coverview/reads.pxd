cdef extern from "htslib/sam.h":
    ctypedef struct bam1_t


cdef class ReadArray:
    cdef bam1_t** reads
    cdef int __size
    cdef int __capacity
    cdef int __longest_read
    cdef void append(self, bam1_t* read)
    cdef void setWindowPointers(self, int start, int end, bam1_t*** window_start, bam1_t*** window_end)
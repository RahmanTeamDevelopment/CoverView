
from pysam.libchtslib cimport bam1_t


cdef class ReadArray:
    cdef bam1_t** reads
    cdef int __size
    cdef int __capacity
    cdef int __longest_read
    cdef void append(self, bam1_t* read)
    cdef void set_pointers_to_start_and_end_of_interval(self, int start, int end, bam1_t*** window_start, bam1_t*** window_end)
    cdef int count_reads_in_interval(self, int start_pos, int end_pos)
"""
Utility classes and functions for efficient processing of read data
"""

from libc.stdint cimport int32_t
from pysam.libcalignmentfile cimport bam1_t

cdef extern from "sam.h":
    bam1_t *bam_dup1(const bam1_t *bsrc)
    void bam_destroy1(bam1_t *b)
    int32_t bam_endpos(const bam1_t *b)

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    int abs(int)

cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  int strncmp(char *s1,char *s2,size_t len)
  char *strncpy(char *dest,char *src, size_t len)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )


cdef int bisectReadsLeft(bam1_t** reads, int testPos, int nReads):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = nReads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if reads[mid].core.pos < testPos:
            low = mid + 1
        else:
            high = mid

    return low


cdef int bisectReadsRight(bam1_t** reads, int testPos, int nReads):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = nReads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if testPos < reads[mid].core.pos:
            high = mid
        else:
            low = mid + 1

    return low


cdef class ReadArray:
    """
    Simple structure to wrap a raw C array, with some bounds
    checking.
    """
    def __init__(self, int size):
        """
        Allocate an array of size 'size', with initial values
        'init'.
        """
        self.reads = <bam1_t**>(malloc(size*sizeof(bam1_t*)))
        assert self.reads != NULL, "Could not allocate memory for ReadArray"
        
        self.__size = 0 # We don't put anything in here yet, just allocate memory
        self.__capacity = size
        self.__longest_read = 0

        cdef int index = 0

        for index from 0 <= index < size:
            self.reads[index] = NULL
    
    def __dealloc__(self):
        """
        Free memory
        """
        cdef int index = 0

        if self.reads != NULL:
            for index from 0 <= index < self.__size:
                if self.reads[index] != NULL:
                    bam_destroy1(self.reads[index])
                    self.reads[index] = NULL
            
            free(self.reads)

    cdef void append(self, bam1_t* read):
        """
        Append a new read to the array, re-allocating if necessary.
        """
        cdef bam1_t** temp = NULL

        if self.__size == self.__capacity:
            temp = <bam1_t**>(realloc(self.reads, 2*sizeof(bam1_t*)*self.__capacity))

            if temp == NULL:
                raise StandardError, "Could not re-allocate ReadArray"
            else:
                self.reads = temp
                self.__capacity *= 2

        self.reads[self.__size] = bam_dup1(read)
        self.__size += 1

        cdef int read_length = bam_endpos(read) - read.core.pos

        if read_length > self.__longest_read:
            self.__longest_read = read_length

    cdef void setWindowPointers(self, int start, int end, bam1_t** window_start, bam1_t** window_end):
        """
        Set the pointers 'windowStart' and 'windowEnd' to
        point to the relevant first and last +1 reads as specified by
        the co-ordinates.
        """
        cdef int firstOverlapStart = -1
        cdef int startPosOfReads = -1
        cdef int endPosOfReads = -1

        if self.__size == 0:
            window_start = self.reads
            window_end = self.reads
        else:
            firstOverlapStart = max(1, start - self.__longest_read)
            startPosOfReads = bisectReadsLeft(self.reads, firstOverlapStart, self.__size)
            endPosOfReads = bisectReadsLeft(self.reads, end, self.__size)

            while startPosOfReads < self.__size and bam_endpos(self.reads[startPosOfReads]) <= start:
                startPosOfReads += 1

            window_start = self.reads + startPosOfReads
            window_end = min(self.reads + endPosOfReads, self.reads + self.__size)

            assert startPosOfReads <= endPosOfReads

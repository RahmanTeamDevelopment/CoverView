

cdef class QualityHistogram:
    cdef int n_data_points
    cdef int* data
    cdef float compute_fraction_below_threshold(self, int threshold)
    cdef float compute_median(self)
    cdef add_data(self, int quality_score)
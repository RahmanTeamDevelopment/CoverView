

cdef class QualityHistogramArray:
    cdef int* n_data_points
    cdef int** data
    cdef int num_hists
    cdef float compute_fraction_below_threshold(self, int index, int threshold)
    cdef float compute_median(self, int index)
    cdef void add_data(self, int index, int quality_score)
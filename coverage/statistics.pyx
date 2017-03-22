"""
Simple statistics utilities for dealing with coverage and read quality
data
"""
from __future__ import division


cdef extern from "stdlib.h":
    void free(void *)
    void* malloc(size_t)
    void* calloc(size_t,size_t)


cdef class QualityHistogramArray:
    """
    Stores an array of histograms of quality scores and computes summary stats
    (mean, median etc) on that histogram. Using histograms for this rather
    than storing the raw quality values is an optimisation for both storage and
    run-time. We store 101s elements per histogram rather than n_bases, and the median
    can be computed in O(N) rather than O(N*logN).
    """
    def __init__(self, int num_hists):
        cdef int size = 101
        self.data = <int**>(malloc(num_hists*sizeof(int*)))
        self.n_data_points = <int*>(calloc(num_hists, sizeof(int)))

        cdef int i = 0

        for i from 0 <= i < num_hists:
            self.data[i] = <int*>(calloc(size, sizeof(int)))

    def __dealloc__(self):
        free(self.data)

    cdef void add_data(self, int index, int quality_score):
        self.n_data_points[index] += 1
        self.data[index][quality_score] += 1

    cdef float compute_fraction_below_threshold(self, int index, int threshold):

        cdef int total = 0
        cdef int i = 0

        if self.n_data_points[index] == 0:
            return float('NaN')
        else:
            for i from 0 <= i < threshold:
                total += self.data[index][i]

            return total / self.n_data_points[index]

    cdef float compute_median(self, int index):

        cdef int total = 0
        cdef int current_bin = 0
        cdef int last_non_empty_bin = -1

        if self.n_data_points[index] == 0:
            return float('NaN')

        if self.n_data_points[index] % 2 == 0:

            while True:
                total += self.data[index][current_bin]

                if total > self.n_data_points[index] // 2:
                    if last_non_empty_bin == -1:
                        return current_bin
                    else:
                        return (current_bin + last_non_empty_bin) / 2.0

                if self.data[index][current_bin] > 0:
                    last_non_empty_bin = current_bin

                current_bin += 1
        else:
            while True:
                total += self.data[index][current_bin]
                if total > self.n_data_points[index] // 2:
                    return current_bin

                current_bin += 1

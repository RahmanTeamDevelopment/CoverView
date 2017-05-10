# cython: profile=False

"""
Simple statistics utilities for dealing with coverage and read quality
data.

N.B. Profiling is disabled for this module. These are all small and highly-optimised 
functions. Profiling results are misleading because the overhead of profiling is significantly
greater than the cost of the function call
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
        self.num_hists = num_hists
        self.data = <int**>(malloc(num_hists*sizeof(int*)))
        self.n_data_points = <int*>(calloc(num_hists, sizeof(int)))

        cdef int i = 0

        for i from 0 <= i < num_hists:
            self.data[i] = <int*>(calloc(size, sizeof(int)))

    def __dealloc__(self):
        cdef int i = 0
        for i from 0 < i < self.num_hists:
            free(self.data[i])

        free(self.data)
        free(self.n_data_points)

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

            return <float>(total) / <float>(self.n_data_points[index])

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


class pyQualityHistogramArray(object):
    """
    Wrapper for the above class. This exists to a) allow us to test the QualityHistogramArray class
    using the python unit-test module and b) expose the functionality of the QualityHistogramArray
    class to pure Python code.
    """
    def __init__(self, int num_hists):
        self._hist_array = QualityHistogramArray(num_hists)

    def add_data(self, int index, int quality_score):
        cdef QualityHistogramArray hist_array = self._hist_array
        return hist_array.add_data(index, quality_score)

    def compute_fraction_below_threshold(self, int index, int threshold):
        cdef QualityHistogramArray hist_array = self._hist_array
        return hist_array.compute_fraction_below_threshold(index, threshold)

    def compute_median(self, int index):
        cdef QualityHistogramArray hist_array = self._hist_array
        return hist_array.compute_median(index)


class pyQualityHistogram(object):
    """
    Simplified wrapper for the QualityHistogramArray class, which specialises the interface for
    the case when there is only a single histogram
    """
    def __init__(self):
        self._hist_array = QualityHistogramArray(1)

    def add_data(self, int quality_score):
        cdef QualityHistogramArray hist_array = self._hist_array
        hist_array.add_data(0, quality_score)

    def compute_fraction_below_threshold(self, int threshold):
        cdef QualityHistogramArray hist_array = self._hist_array
        return hist_array.compute_fraction_below_threshold(0, threshold)

    def compute_median(self):
        cdef QualityHistogramArray hist_array = self._hist_array
        return hist_array.compute_median(0)


def median(list x):
    """
    Calculate and return the median value of an input list. The median is the central value.    
    """
    list_length = len(x)

    if list_length == 0:
        return float('NaN')

    sorted_list = sorted(x)

    if list_length % 2 == 0:
        lower_index = (list_length // 2) - 1
        upper_index = list_length // 2
        return 0.5 * (
            sorted_list[lower_index] + sorted_list[upper_index]
        )
    else:
        return sorted_list[list_length // 2]
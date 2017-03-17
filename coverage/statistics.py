"""
Simple statistics utilities for dealing with coverage and read quality
data
"""
from __future__ import division


class QualityHistogram(object):
    """
    Stores a histogram of quality scores and computes summary stats
    (mean, median etc) on that histogram. Using a histogram for this rather
    than storing the raw quality values is an optimisation for both storage and
    run-time. We store 101s elements per histogram rather than n_bases, and the median
    can be computed in O(N) rather than O(N*logN).
    """
    def __init__(self):
        self.data = [0]*101
        self.n_data_points = 0

    def add_data(self, quality_score):
        self.n_data_points += 1
        self.data[quality_score] += 1

    def compute_fraction_below_threshold(self, threshold):

        if self.n_data_points == 0:
            return float('NaN')
        else:
            total = 0

            for i in xrange(threshold):
                total += self.data[i]

            return total / self.n_data_points

    def compute_median(self):

        if self.n_data_points == 0:
            return float('NaN')

        total = 0
        current_bin = 0

        if self.n_data_points % 2 == 0:
            last_non_empty_bin = None

            while True:
                total += self.data[current_bin]

                if total > self.n_data_points // 2:
                    return (current_bin + last_non_empty_bin) / 2.0

                if self.data[current_bin] > 0:
                    last_non_empty_bin = current_bin

                current_bin += 1
        else:
            while True:
                total += self.data[current_bin]
                if total > self.n_data_points // 2:
                    return current_bin

                current_bin += 1

import coverview
import math
import unittest


class TestMedianCalculation(unittest.TestCase):

    def test_empty_histogram_has_median_of_None(self):
        hist = coverview.statistics.pyQualityHistogram()
        assert math.isnan(hist.compute_median())

    def test_median_of_single_value_is_that_value(self):
        hist = coverview.statistics.pyQualityHistogram()
        hist.add_data(10)
        assert hist.compute_median() == 10

    def test_median_of_two_values_is_the_mean_of_those_values(self):
        hist = coverview.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(20)
        assert hist.compute_median() == 15

    def test_median_of_two_identical_values_is_the_value(self):
        hist = coverview.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(10)
        assert hist.compute_median() == 10

    def test_median_of_three_values_is_the_middle_value(self):
        hist = coverview.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(20)
        hist.add_data(30)
        assert hist.compute_median() == 20

    def test_median_of_0_to_99_is_49_and_a_half(self):
        hist = coverview.statistics.pyQualityHistogram()
        for i in xrange(0, 100):
            hist.add_data(i)
        assert hist.compute_median() == 49.5

    def test_median_of_0_to_100_is_50(self):
        hist = coverview.statistics.pyQualityHistogram()
        for i in xrange(0, 101):
            hist.add_data(i)
        assert hist.compute_median() == 50.0


class TestFractionBelowThresholdCalculation(object):

    def test_fraction_is_nan_for_empty_histogram(self):
        hist = coverview.statistics.pyQualityHistogram()
        assert math.isnan(hist.compute_fraction_below_threshold(100))

    def test_fraction_is_zero_for_single_value_above_threshold(self):
        hist = coverview.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold + 1)
        assert hist.compute_fraction_below_threshold(threshold) == 0.0

    def test_fraction_is_zero_for_single_value_at_threshold(self):
        hist = coverview.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold)
        assert hist.compute_fraction_below_threshold(threshold) == 0.0

    def test_fraction_is_one_for_single_value_below_threshold(self):
        hist = coverview.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold + 1)
        assert hist.compute_fraction_below_threshold(threshold) == 0.0

    def test_fraction_is_one_half_for_two_values_one_above_and_one_below_threshold(self):
        hist = coverview.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold + 1)
        hist.add_data(threshold - 1)
        assert hist.compute_fraction_below_threshold(threshold) == 0.5


if __name__ == "__main__":
    unittest.main()



import coverview_.statistics
import math
import unittest


class TestSimpleMedianCalculation(unittest.TestCase):

    def test_empty_list_has_median_of_nan(selfs):
        x = []
        assert math.isnan(coverview_.statistics.median(x))

    def test_median_of_single_value_is_that_value(self):
        x = [5]
        assert coverview_.statistics.median(x) == 5

    def test_median_of_two_values_is_the_mean_of_those_values(self):
        x = [10, 20]
        assert coverview_.statistics.median(x) == 15

    def test_median_of_two_identical_values_is_the_value(self):
        x = [10, 10]
        assert coverview_.statistics.median(x) == 10

    def test_median_of_three_values_is_the_middle_value(self):
        x = [10, 20, 30]
        assert coverview_.statistics.median(x) == 20

    def test_median_of_0_to_99_is_49_and_a_half(self):
        x = range(0, 100)
        assert coverview_.statistics.median(x) == 49.5

    def test_median_of_0_to_100_is_50(self):
        x = range(0, 101)
        assert coverview_.statistics.median(x) == 50.0


class TestQualityHistogramMedianCalculation(unittest.TestCase):

    def test_empty_histogram_has_median_of_nan(self):
        hist = coverview_.statistics.pyQualityHistogram()
        assert math.isnan(hist.compute_median())

    def test_median_of_single_value_is_that_value(self):
        hist = coverview_.statistics.pyQualityHistogram()
        hist.add_data(10)
        assert hist.compute_median() == 10

    def test_median_of_two_values_is_the_mean_of_those_values(self):
        hist = coverview_.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(20)
        assert hist.compute_median() == 15

    def test_median_of_two_identical_values_is_the_value(self):
        hist = coverview_.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(10)
        assert hist.compute_median() == 10

    def test_median_of_three_values_is_the_middle_value(self):
        hist = coverview_.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(20)
        hist.add_data(30)
        assert hist.compute_median() == 20

    def test_median_of_four_values_with_one_low_and_three_high_is_high_value(self):
        hist = coverview_.statistics.pyQualityHistogram()
        hist.add_data(10)
        hist.add_data(20)
        hist.add_data(20)
        hist.add_data(20)
        assert hist.compute_median() == 20

    def test_median_of_0_to_99_is_49_and_a_half(self):
        hist = coverview_.statistics.pyQualityHistogram()
        for i in range(0, 100):
            hist.add_data(i)
        assert hist.compute_median() == 49.5

    def test_median_of_0_to_100_is_50(self):
        hist = coverview_.statistics.pyQualityHistogram()
        for i in range(0, 101):
            hist.add_data(i)
        assert hist.compute_median() == 50.0


class TestQualityHistogramFractionBelowThresholdCalculation(object):

    def test_fraction_is_nan_for_empty_histogram(self):
        hist = coverview_.statistics.pyQualityHistogram()
        assert math.isnan(hist.compute_fraction_below_threshold(100))

    def test_fraction_is_zero_for_single_value_above_threshold(self):
        hist = coverview_.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold + 1)
        assert hist.compute_fraction_below_threshold(threshold) == 0.0

    def test_fraction_is_zero_for_single_value_at_threshold(self):
        hist = coverview_.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold)
        assert hist.compute_fraction_below_threshold(threshold) == 0.0

    def test_fraction_is_one_for_single_value_below_threshold(self):
        hist = coverview_.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold + 1)
        assert hist.compute_fraction_below_threshold(threshold) == 0.0

    def test_fraction_is_one_half_for_two_values_one_above_and_one_below_threshold(self):
        hist = coverview_.statistics.pyQualityHistogram()
        threshold = 10
        hist.add_data(threshold + 1)
        hist.add_data(threshold - 1)
        assert hist.compute_fraction_below_threshold(threshold) == 0.5


if __name__ == "__main__":
    unittest.main()

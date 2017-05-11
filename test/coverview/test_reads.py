import coverview
import math
import unittest


class TestReadArray(unittest.TestCase):

    def test_empty_list_has_median_of_nan(selfs):
        x = []
        assert math.isnan(coverview.statistics.median(x))

    def test_median_of_single_value_is_that_value(self):
        x = [5]
        assert coverview.statistics.median(x) == 5

    def test_median_of_two_values_is_the_mean_of_those_values(self):
        x = [10, 20]
        assert coverview.statistics.median(x) == 15

    def test_median_of_two_identical_values_is_the_value(self):
        x = [10, 10]
        assert coverview.statistics.median(x) == 10

    def test_median_of_three_values_is_the_middle_value(self):
        x = [10, 20, 30]
        assert coverview.statistics.median(x) == 20

    def test_median_of_0_to_99_is_49_and_a_half(self):
        x = range(0, 100)
        assert coverview.statistics.median(x) == 49.5

    def test_median_of_0_to_100_is_50(self):
        x = range(0, 101)
        assert coverview.statistics.median(x) == 50.0


if __name__ == "__main__":
    unittest.main()
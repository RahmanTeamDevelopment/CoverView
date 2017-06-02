import coverview.utils
import math
import unittest


class TestMinOrNaN(unittest.TestCase):
    def test_returns_nan_if_empty_list_is_passed(self):
        assert math.isnan(
            coverview.utils.min_or_nan([])
        )

    def test_returns_nan_if_nan_is_passed(self):
        assert math.isnan(
            coverview.utils.min_or_nan([float('NaN')])
        )

    def test_returns_value_if_single_value_is_passed(self):
        assert coverview.utils.min_or_nan([10]) == 10

    def test_returns_value_if_single_value_and_nan_are_passed(self):
        assert coverview.utils.min_or_nan([10, float('NaN')]) == 10

    def test_returns_smaller_value_if_two_values_are_passed(self):
        assert coverview.utils.min_or_nan([10, 20]) == 10


class TestMaxOrNaN(unittest.TestCase):
    def test_returns_nan_if_empty_list_is_passed(self):
        assert math.isnan(
            coverview.utils.max_or_nan([])
        )

    def test_returns_nan_if_nan_is_passed(self):
        assert math.isnan(
            coverview.utils.max_or_nan([float('NaN')])
        )

    def test_returns_value_if_single_value_is_passed(self):
        assert coverview.utils.max_or_nan([10]) == 10

    def test_returns_value_if_single_value_and_nan_are_passed(self):
        assert coverview.utils.max_or_nan([10, float('NaN')]) == 10

    def test_returns_larger_value_if_two_values_are_passed(self):
        assert coverview.utils.max_or_nan([10, 20]) == 20


class TestGetNamesOfTargetRegions(unittest.TestCase):
    def test_raises_exception_if_bed_file_is_empty(self):
        self.assertRaises(
            StandardError,
            coverview.utils.get_names_of_target_regions,
            []
        )

    def test_with_single_region(self):
        names, num_unique_names = coverview.utils.get_names_of_target_regions([
            "1\t32\t32\tRegion_1"
        ])

        assert len(names) == 1
        assert num_unique_names == 1
        assert names[0] == "Region_1"

    def test_with_two_regions_with_different_names(self):
        names, num_unique_names = coverview.utils.get_names_of_target_regions([
            "1\t32\t32\tRegion_1",
            "1\t1000\t1100\tRegion_2"
        ])

        assert len(names) == 2
        assert num_unique_names == 2
        assert names[0] == "Region_1"
        assert names[1] == "Region_2"

    def test_with_two_regions_with_the_same_name(self):
        names, num_unique_names = coverview.utils.get_names_of_target_regions([
            "1\t32\t32\tRegion_1",
            "1\t1000\t1100\tRegion_1"
        ])

        assert len(names) == 2
        assert num_unique_names == 1
        assert names[0] == "Region_1_1"
        assert names[1] == "Region_1_2"


class Tes
        # def get_clusters_of_regions_from_bed_file(bed_file_name, size_limit=100000):


if __name__ == "__main__":
    unittest.main()
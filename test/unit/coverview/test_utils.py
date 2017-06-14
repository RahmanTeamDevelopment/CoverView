import copy
import coverview.utils
import math
import unittest


class TestUniquifyRegionNames(unittest.TestCase):
    def test_returns_empty_list_when_given_empty_list(self):
        assert coverview.utils.uniquify_region_names([]) == []

    def test_returns_list_of_one_interval_when_given_list_of_one_interval(self):
        interval = coverview.utils.GenomicInterval("1", 100, 200, "REGION")
        assert coverview.utils.uniquify_region_names([interval]) == [interval]

    def test_returns_identical_copy_of_inputs_when_inputs_have_different_names(self):
        interval_1 = coverview.utils.GenomicInterval("1", 100, 200, "REGION")
        interval_2 = coverview.utils.GenomicInterval("1", 500, 600, "OTHER_REGION")
        inputs = [interval_1, interval_2]
        assert coverview.utils.uniquify_region_names(inputs) == inputs

    def test_returns_copy_with_unique_names_if_identical_names_in_input(self):
        interval_1 = coverview.utils.GenomicInterval("1", 100, 200, "REGION")
        interval_2 = coverview.utils.GenomicInterval("1", 500, 600, "REGION")
        inputs = [interval_1, interval_2]

        expected_outputs = [
            coverview.utils.GenomicInterval("1", 100, 200, "REGION_1"),
            coverview.utils.GenomicInterval("1", 500, 600, "REGION_2"),
        ]

        assert coverview.utils.uniquify_region_names(inputs) == expected_outputs

    def test_returns_sorted_copy_with_unique_names_if_identical_names_in_input_and_inputs_unsorted(self):
        interval_1 = coverview.utils.GenomicInterval("1", 500, 600, "REGION")
        interval_2 = coverview.utils.GenomicInterval("1", 100, 200, "REGION")
        inputs = [interval_1, interval_2]

        expected_outputs = [
            coverview.utils.GenomicInterval("1", 100, 200, "REGION_1"),
            coverview.utils.GenomicInterval("1", 500, 600, "REGION_2"),
        ]

        assert coverview.utils.uniquify_region_names(inputs) == expected_outputs


class TestGenpmicInterval(unittest.TestCase):
    def test_raises_assertion_if_constructed_with_negative_position(self):
        self.assertRaises(
            AssertionError,
            coverview.utils.GenomicInterval,
            chromosome="1",
            start_pos=-1,
            end_pos=100,
            name=None
        )

        self.assertRaises(
            AssertionError,
            coverview.utils.GenomicInterval,
            chromosome="1",
            start_pos=100,
            end_pos=-1,
            name=None
        )

    def test_interval_size_is_end_minus_start(self):
        interval_1 = coverview.utils.GenomicInterval(None, 100, 200)
        interval_2 = coverview.utils.GenomicInterval(None, 0, 0)
        assert interval_1.size() == 100
        assert interval_2.size() == 0

    def test_intervals_compare_equal_if_chromosomes_and_coordinates_are_identical(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)
        interval_2 = coverview.utils.GenomicInterval("1", 0, 100)
        assert interval_1 == interval_2

    def test_intervals_do_not_compare_equal_if_chromosomes_are_different(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)
        interval_2 = coverview.utils.GenomicInterval("2", 0, 100)
        assert interval_1 != interval_2

    def test_intervals_do_not_compare_equal_if_start_positions_are_different(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)
        interval_2 = coverview.utils.GenomicInterval("1", 10, 100)
        assert interval_1 != interval_2

    def test_intervals_do_not_compare_equal_if_end_positions_are_different(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)
        interval_2 = coverview.utils.GenomicInterval("1", 0, 150)
        assert interval_1 != interval_2

    def test_overlap_of_identical_intervals_is_same_interval_again(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)
        interval_2 = copy.deepcopy(interval_1)
        assert interval_1.overlap(interval_2) == interval_1

    def test_overlap_of_intervals_on_different_chromosomes_is_zero(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)
        interval_2 = coverview.utils.GenomicInterval("2", 0, 100)
        assert interval_1.overlap(interval_2).size() == 0

    def test_overlap_of_non_overlapping_intervals_on_same_chromosome_is_zero(self):
        interval_1 = coverview.utils.GenomicInterval("1", 0, 100)

        interval_2 = coverview.utils.GenomicInterval(
            interval_1.chromosome,
            interval_1.end_pos + 1,
            interval_1.end_pos + 2
        )

        assert interval_1.overlap(interval_2).size() == 0

    def test_overlap_of_small_interval_entirely_contained_in_larger_interval_is_small_interval(self):
        small_interval = coverview.utils.GenomicInterval("1", 100, 200)

        larger_interval = coverview.utils.GenomicInterval(
            small_interval.chromosome,
            small_interval.start_pos - 100,
            small_interval.end_pos + 100
        )

        assert small_interval.overlap(larger_interval) == small_interval


class TestBEDParser(unittest.TestCase):
    def test_raises_standard_error_if_line_in_file_has_less_than_three_columns(self):
        parser = coverview.utils.BedFileParser(iter([
            "chr1\t100\n"
        ]))

        self.assertRaises(
            StandardError,
            parser.next
        )

    def test_returns_genomic_interval_from_well_formatted_line(self):
        parser = coverview.utils.BedFileParser(iter([
            "chr1\t100\t200\tREGION_1\n",
            "chr1\t300\t400\tREGION_2\n",
        ]))

        assert parser.next() == coverview.utils.GenomicInterval(
            "chr1",
            100,
            200,
            "REGION_1"
        )

        assert parser.next() == coverview.utils.GenomicInterval(
            "chr1",
            300,
            400,
            "REGION_2"
        )


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
    def test_raises_exception_if_there_are_no_regions(self):
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


class TestGetClustersOfRegionsFromBEDFile(unittest.TestCase):
    def test_raises_exception_if_there_are_no_regions(self):

        # Generator functiond won't raise an exception until we start iterating
        # through the results
        clusters = coverview.utils.get_clusters_of_regions_from_bed_file(
            []
        )

        self.assertRaises(
            StandardError,
            clusters.next
        )


if __name__ == "__main__":
    unittest.main()
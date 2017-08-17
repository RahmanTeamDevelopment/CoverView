import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverViewRegionsOutput(unittest.TestCase):

    def test_regions_output_with_zero_coverage_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"

            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132

            assert regions_output['Region_1']["RC"] == 0
            assert regions_output['Region_1']["RC+"] == 0
            assert regions_output['Region_1']["RC-"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 0
            assert regions_output['Region_1']["MINCOV"] == 0
            assert regions_output['Region_1']["MEDQCOV"] == 0
            assert regions_output['Region_1']["MINQCOV"] == 0

            assert regions_output['Region_1']["MAXFLMQ"] == "."
            assert regions_output['Region_1']["MAXFLBQ"] == "."

    def test_regions_output_one_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_1']["RC"] == 1
            assert regions_output['Region_1']["RC+"] == 1
            assert regions_output['Region_1']["RC-"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 1
            assert regions_output['Region_1']["MINCOV"] == 1
            assert regions_output['Region_1']["MEDQCOV"] == 1
            assert regions_output['Region_1']["MINQCOV"] == 1
            assert regions_output['Region_1']["MAXFLMQ"] == 0.0
            assert regions_output['Region_1']["MAXFLBQ"] == 0.0

    def test_regions_output_two_reads_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 2))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_1']["RC"] == 2
            assert regions_output['Region_1']["RC+"] == 2
            assert regions_output['Region_1']["RC-"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 2
            assert regions_output['Region_1']["MINCOV"] == 2
            assert regions_output['Region_1']["MEDQCOV"] == 2
            assert regions_output['Region_1']["MINQCOV"] == 2
            assert regions_output['Region_1']["MAXFLMQ"] == 0.0
            assert regions_output['Region_1']["MAXFLBQ"] == 0.0

    def test_regions_output_many_reads_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 100))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_1']["RC"] == 100
            assert regions_output['Region_1']["RC+"] == 100
            assert regions_output['Region_1']["RC-"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 100
            assert regions_output['Region_1']["MINCOV"] == 100
            assert regions_output['Region_1']["MEDQCOV"] == 100
            assert regions_output['Region_1']["MINQCOV"] == 100
            assert regions_output['Region_1']["MAXFLMQ"] == 0.0
            assert regions_output['Region_1']["MAXFLBQ"] == 0.0

    def test_regions_output_with_two_regions(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_reads(("1", 1000, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_region(("1", 1000, 1100, "Region_2"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert "Region_2" in regions_output

            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_1']["RC"] == 1
            assert regions_output['Region_1']["RC+"] == 1
            assert regions_output['Region_1']["RC-"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 1
            assert regions_output['Region_1']["MINCOV"] == 1
            assert regions_output['Region_1']["MEDQCOV"] == 1
            assert regions_output['Region_1']["MINQCOV"] == 1
            assert regions_output['Region_1']["MAXFLMQ"] == 0.0
            assert regions_output['Region_1']["MAXFLBQ"] == 0.0

            assert regions_output['Region_2']["Chromosome"] == "1"
            assert regions_output['Region_2']["Start_position"] == 1000
            assert regions_output['Region_2']["End_position"] == 1100
            assert regions_output['Region_2']["RC"] == 1
            assert regions_output['Region_2']["RC+"] == 1
            assert regions_output['Region_2']["RC-"] == 0
            assert regions_output['Region_2']["MEDCOV"] == 1
            assert regions_output['Region_2']["MINCOV"] == 1
            assert regions_output['Region_2']["MEDQCOV"] == 1
            assert regions_output['Region_2']["MINQCOV"] == 1
            assert regions_output['Region_2']["MAXFLMQ"] == 0.0
            assert regions_output['Region_2']["MAXFLBQ"] == 0.0

    def test_regions_output_with_two_adjacent_regions_of_different_coverage(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_reads(("1", 132, 100, 10))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_region(("1", 132, 232, "Region_2"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert "Region_2" in regions_output

            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_1']["RC"] == 1
            assert regions_output['Region_1']["RC+"] == 1
            assert regions_output['Region_1']["RC-"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 1
            assert regions_output['Region_1']["MINCOV"] == 1
            assert regions_output['Region_1']["MEDQCOV"] == 1
            assert regions_output['Region_1']["MINQCOV"] == 1
            assert regions_output['Region_1']["MAXFLMQ"] == 0.0
            assert regions_output['Region_1']["MAXFLBQ"] == 0.0

            assert regions_output['Region_2']["Chromosome"] == "1"
            assert regions_output['Region_2']["Start_position"] == 132
            assert regions_output['Region_2']["End_position"] == 232
            assert regions_output['Region_2']["RC"] == 10
            assert regions_output['Region_2']["RC+"] == 10
            assert regions_output['Region_2']["RC-"] == 0
            assert regions_output['Region_2']["MEDCOV"] == 10
            assert regions_output['Region_2']["MINCOV"] == 10
            assert regions_output['Region_2']["MEDQCOV"] == 10
            assert regions_output['Region_2']["MINQCOV"] == 10
            assert regions_output['Region_2']["MAXFLMQ"] == 0.0
            assert regions_output['Region_2']["MAXFLBQ"] == 0.0

    def test_that_regions_with_same_name_have_1_based_index_appended_in_output_file(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_reads(("1", 132, 100, 10))
            runner.add_region(("1", 32, 132, "Region"))
            runner.add_region(("1", 132, 232, "Region"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert "Region_2" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 32
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_2']["Chromosome"] == "1"
            assert regions_output['Region_2']["Start_position"] == 132
            assert regions_output['Region_2']["End_position"] == 232

    def test_that_when_input_regions_are_unsorted_output_regions_are_sorted(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 100, 100, 1))
            runner.add_reads(("1", 500, 100, 10))
            runner.add_region(("1", 500, 600, "Region"))
            runner.add_region(("1", 100, 200, "Region"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert "Region_2" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 100
            assert regions_output['Region_1']["End_position"] == 200
            assert regions_output['Region_2']["Chromosome"] == "1"
            assert regions_output['Region_2']["Start_position"] == 500
            assert regions_output['Region_2']["End_position"] == 600


if __name__ == "__main__":
    unittest.main()

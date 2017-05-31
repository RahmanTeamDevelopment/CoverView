import math
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

            regions_output = testutils.output_checkers.load_coveriew_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"

            # Output positions are 1-indexed and intervals are closed
            assert regions_output['Region_1']["Start_position"] == 33
            assert regions_output['Region_1']["End_position"] == 132

            # No data so all coverages are 0
            assert regions_output['Region_1']["RC"] == 0
            assert regions_output['Region_1']["MEDCOV"] == 0
            assert regions_output['Region_1']["MINCOV"] == 0
            assert regions_output['Region_1']["MEDQCOV"] == 0
            assert regions_output['Region_1']["MINQCOV"] == 0

            # Fractions are NaN with 0 reads
            assert math.isnan(regions_output['Region_1']["MAXFLMQ"])
            assert math.isnan(regions_output['Region_1']["MAXFLBQ"])

    def test_regions_output_one_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coveriew_regions_output(
                "output_regions.txt"
            )

            assert "Region_1" in regions_output
            assert regions_output['Region_1']["Chromosome"] == "1"
            assert regions_output['Region_1']["Start_position"] == 33
            assert regions_output['Region_1']["End_position"] == 132
            assert regions_output['Region_1']["RC"] == 1
            assert regions_output['Region_1']["MEDCOV"] == 1
            assert regions_output['Region_1']["MINCOV"] == 1
            assert regions_output['Region_1']["MEDQCOV"] == 1
            assert regions_output['Region_1']["MINQCOV"] == 1
            assert regions_output['Region_1']["MAXFLMQ"] == 0.0
            assert regions_output['Region_1']["MAXFLBQ"] == 0.0


if __name__ == "__main__":
    unittest.main()
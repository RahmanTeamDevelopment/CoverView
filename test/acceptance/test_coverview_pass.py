import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverViewWithPassCriteria(unittest.TestCase):

    def test_regions_with_minqcov_below_min_threshold_are_marked_as_fail(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },

                "pass": {
                    "MINQCOV_MIN": 1
                },

                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert regions_output['Region_1']['Pass_or_fail'] == "FAIL"

    def test_regions_with_minqcov_above_min_threshold_are_marked_as_pass(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 2))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },

                "pass": {
                    "MINQCOV_MIN": 1
                },

                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert regions_output['Region_1']['Pass_or_fail'] == "PASS"

    def test_regions_with_maxqcov_below_max_threshold_are_marked_as_pass(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },

                "pass": {
                    "MINQCOV_MAX": 1
                },

                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert regions_output['Region_1']['Pass_or_fail'] == "PASS"

    def test_regions_with_maxqcov_above_max_threshold_are_marked_as_fail(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 2))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },

                "pass": {
                    "MINQCOV_MAX": 1
                },

                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coverview_regions_output(
                "output_regions.txt"
            )

            assert regions_output['Region_1']['Pass_or_fail'] == "FAIL"

    def test_coverview_raises_exception_when_invalid_pass_fail_config_is_used(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },

                "fail": {
                    "MIN_MINQCOV": 50,
                    "MAX_MAXFLMQ": 0.05,
                    "MAX_MAXFLBQ": 0.15
                },
            })

            self.assertRaises(
                StandardError,
                runner.run_coverview_and_get_exit_code
            )


if __name__ == "__main__":
    unittest.main()

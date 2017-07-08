import bamgen.bamgen
import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverviewWithQualityThresholds(unittest.TestCase):

    def test_bases_with_base_quality_less_than_threshold_not_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_bq": bamgen.bamgen.HIGH_QUALITY + 1
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 0

    def test_bases_with_base_quality_equal_to_threshold_are_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_bq": bamgen.bamgen.HIGH_QUALITY
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 1

    def test_bases_with_base_quality_greater_than_threshold_are_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_bq": bamgen.bamgen.HIGH_QUALITY - 1
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 1

    def test_reads_with_mapping_quality_less_than_threshold_not_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_mq": bamgen.bamgen.HIGH_QUALITY + 1
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 0

    def test_reads_with_mapping_quality_equal_to_threshold_are_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_mq": bamgen.bamgen.HIGH_QUALITY
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 1

    def test_reads_with_mapping_quality_greater_than_threshold_are_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_mq": bamgen.bamgen.HIGH_QUALITY - 1
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 1

    def test_reads_with_mapping_and_base_quality_above_threshold_are_counted_in_qcov(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },
                "low_bq": bamgen.bamgen.HIGH_QUALITY - 1,
                "low_mq": bamgen.bamgen.HIGH_QUALITY - 1
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert profile_output['Region_1']["1:36"]['COV'] == 1
            assert profile_output['Region_1']["1:36"]['QCOV'] == 1


if __name__ == "__main__":
    unittest.main()

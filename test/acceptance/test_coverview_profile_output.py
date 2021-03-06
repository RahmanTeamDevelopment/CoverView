import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverViewProfileOutput(unittest.TestCase):

    def test_with_zero_coverage_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({"outputs": {"profiles": True}})
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" in profile_output
            assert "1:40" not in profile_output['Region_1']

            for position in [35, 36, 37, 38, 39]:
                chrom_pos = "1:{}".format(position)
                assert profile_output['Region_1'][chrom_pos]['#Chromosome'] == "1"
                assert profile_output['Region_1'][chrom_pos]['Position'] == position
                assert profile_output['Region_1'][chrom_pos]['COV'] == 0
                assert profile_output['Region_1'][chrom_pos]['QCOV'] == 0
                assert profile_output['Region_1'][chrom_pos]['MEDBQ'] == "."
                assert profile_output['Region_1'][chrom_pos]['FLBQ'] == "."
                assert profile_output['Region_1'][chrom_pos]['MEDMQ'] == "."
                assert profile_output['Region_1'][chrom_pos]['FLMQ'] == "."

    def test_with_one_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({"outputs": {"profiles": True}})
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" in profile_output
            assert "1:40" not in profile_output['Region_1']

            for position in [35, 36, 37, 38, 39]:
                chrom_pos = "1:{}".format(position)
                assert profile_output['Region_1'][chrom_pos]['#Chromosome'] == "1"
                assert profile_output['Region_1'][chrom_pos]['Position'] == position
                assert profile_output['Region_1'][chrom_pos]['COV'] == 1
                assert profile_output['Region_1'][chrom_pos]['QCOV'] == 1
                assert profile_output['Region_1'][chrom_pos]['MEDBQ'] == 60.0
                assert profile_output['Region_1'][chrom_pos]['FLBQ'] == 0.0
                assert profile_output['Region_1'][chrom_pos]['MEDMQ'] == 60.0
                assert profile_output['Region_1'][chrom_pos]['FLMQ'] == 0.0

    def test_with_one_low_mapping_quality_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({"outputs": {"profiles": True}})
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" in profile_output
            assert "1:40" not in profile_output['Region_1']

            for position in [35, 36, 37, 38, 39]:
                chrom_pos = "1:{}".format(position)
                assert profile_output['Region_1'][chrom_pos]['#Chromosome'] == "1"
                assert profile_output['Region_1'][chrom_pos]['Position'] == position
                assert profile_output['Region_1'][chrom_pos]['COV'] == 1
                assert profile_output['Region_1'][chrom_pos]['QCOV'] == 1
                assert profile_output['Region_1'][chrom_pos]['MEDBQ'] == 60.0
                assert profile_output['Region_1'][chrom_pos]['FLBQ'] == 0.0
                assert profile_output['Region_1'][chrom_pos]['MEDMQ'] == 60.0
                assert profile_output['Region_1'][chrom_pos]['FLMQ'] == 0.0


class TestOnlyFailProfilesOptionWithPassCriteria(unittest.TestCase):

    def test_fail_regions_are_in_profile_output_when_only_fail_profiles_is_true(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({
                "outputs": {"profiles": True},
                "only_fail_profiles": True,
                "pass": {
                    "MINQCOV_MIN": 1
                },
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" in profile_output

    def test_pass_regions_are_not_in_profile_output_when_only_fail_profiles_is_true(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 2))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({
                "outputs": {"profiles": True},
                "only_fail_profiles": True,
                "pass": {
                    "MINQCOV_MIN": 1
                },
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" not in profile_output

    def test_fail_regions_are_in_profile_output_when_only_fail_profiles_is_false(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({
                "outputs": {"profiles": True},
                "only_fail_profiles": False,
                "pass": {
                    "MINQCOV_MIN": 1
                },
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" in profile_output

    def test_pass_regions_are_in_profile_output_when_only_fail_profiles_is_false(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 5))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_config_data({
                "outputs": {"profiles": True},
                "only_fail_profiles": False,
                "pass": {
                    "MINQCOV_MIN": 1
                },
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            profile_output = testutils.output_checkers.load_coverview_profile_output(
                "output_profiles.txt"
            )

            assert "Region_1" in profile_output


if __name__ == "__main__":
    unittest.main()

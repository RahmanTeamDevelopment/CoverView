import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverViewRegionsOutput(unittest.TestCase):

    def test_summary_output_with_zero_coverage_in_bam(self):
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

            for position in ['36', '37', '38', '39', '40']:
                chrom_pos = "1:{}".format(position)
                assert profile_output['Region_1'][chrom_pos]['#Chromosome'] == "1"
                assert profile_output['Region_1'][chrom_pos]['Position'] == position
                assert profile_output['Region_1'][chrom_pos]['COV'] == "0"
                assert profile_output['Region_1'][chrom_pos]['QCOV'] == "0"
                assert profile_output['Region_1'][chrom_pos]['MEDBQ'] == "."
                assert profile_output['Region_1'][chrom_pos]['FLBQ'] == "."
                assert profile_output['Region_1'][chrom_pos]['MEDMQ'] == "."
                assert profile_output['Region_1'][chrom_pos]['FLMQ'] == "."

    def test_summary_output_with_one_read_in_bam(self):
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

            for position in ['36', '37', '38', '39', '40']:
                chrom_pos = "1:{}".format(position)
                assert profile_output['Region_1'][chrom_pos]['#Chromosome'] == "1"
                assert profile_output['Region_1'][chrom_pos]['Position'] == position
                assert profile_output['Region_1'][chrom_pos]['COV'] == "1"
                assert profile_output['Region_1'][chrom_pos]['QCOV'] == "1"
                assert profile_output['Region_1'][chrom_pos]['MEDBQ'] == "60.0"
                assert profile_output['Region_1'][chrom_pos]['FLBQ'] == "0.0"
                assert profile_output['Region_1'][chrom_pos]['MEDMQ'] == "60.0"
                assert profile_output['Region_1'][chrom_pos]['FLMQ'] == "0.0"


if __name__ == "__main__":
    unittest.main()
import coverview.main
import os
import testutils.runners
import unittest
import uuid


class TestCoverViewSkippingDuplicateReads(unittest.TestCase):

    def setUp(self):
        self.unique_bam_file_name = str(uuid.uuid4()) + ".bam"
        self.unique_index_file_name = self.unique_bam_file_name + ".bai"
        self.unique_bed_file_name = self.unique_bam_file_name.replace(".bam", ".bed")
        self.unique_config_file_name = self.unique_bam_file_name.replace(".bam", ".json")

    def tearDown(self):
        os.remove(self.unique_bam_file_name)
        os.remove(self.unique_index_file_name)
        os.remove(self.unique_bed_file_name)
        os.remove(self.unique_config_file_name)
        os.remove("output_regions.txt")
        os.remove("output_profiles.txt")
        os.remove("output_summary.txt")

    def test_coverview_runs_and_returns_0_exit_code(self):
        read_sets = [
            ("1", 32, 100, 0)
        ]

        regions = [
            ( "1", 32, 132, "Region_1")
        ]

        config = \
            {
                "outputs": {
                    "profiles": True,
                },

                "pass": {
                    "MINQCOV_MIN": 15
                },

                "count_duplicate_reads": False,
                "direction": True
            }

        testutils.runners.bamgen.make_bam_file(self.unique_bam_file_name, read_sets)
        testutils.runners.make_bed_file(self.unique_bed_file_name, regions)
        testutils.runners.make_config_file(self.unique_config_file_name, config)

        command_line_args = testutils.runners.make_command_line_arguments(
            bam_file_name=self.unique_bam_file_name,
            bed_file_name=self.unique_bed_file_name,
            reference_file_name="__MOCK__",
            config_file_name=self.unique_config_file_name
        )

        assert coverview.main.main(command_line_args) == 0


if __name__ == "__main__":
    unittest.main()
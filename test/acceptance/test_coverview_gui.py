import coverview.main
import os
import shutil
import testutils.coverview
import unittest
import uuid


class TestCoverViewWithGuiOutput(unittest.TestCase):

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
        shutil.rmtree("output_gui")

    def test_coverview_runs_and_returns_0_exit_code(self):
        read_sets = [
            ("1", 32, 100, 0)
        ]

        regions = [
            ( "1", 32, 132, "Region_1")
        ]

        coverview_dir = os.getcwd()
        gui_dir = os.path.join(coverview_dir, "gui")

        config = \
            {
                "outputs": {
                    "gui": True,
                    "gui_output_directory": os.path.join(coverview_dir, "output_gui")
                },

                "gui": {
                    "template_gui_html_file": os.path.join(gui_dir, "gui.html"),
                    "javascript_directory": os.path.join(gui_dir, "lib")
                }
            }

        testutils.coverview.bamgen.make_bam_file(self.unique_bam_file_name, read_sets)
        testutils.coverview.make_bed_file(self.unique_bed_file_name, regions)
        testutils.coverview.make_config_file(self.unique_config_file_name, config)

        command_line_args = testutils.coverview.make_command_line_arguments(
            bam_file_name=self.unique_bam_file_name,
            bed_file_name=self.unique_bed_file_name,
            reference_file_name="__MOCK__",
            config_file_name=self.unique_config_file_name
        )

        assert coverview.main.main(command_line_args) == 0


if __name__ == "__main__":
    unittest.main()
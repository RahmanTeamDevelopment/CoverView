import testutils.output_checkers
import testutils.runners
import unittest


class TestGenerationOfGuiOutputJSON(unittest.TestCase):

    def test_generates_json_file_with_correct_data(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            gui_output_file = "output_gui_file.json"
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_gui_output_file(gui_output_file)
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                    "regions": True
                },
                "pass": {
                    "MINQCOV_MIN": 15
                },
                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0
            json_data = testutils.output_checkers.load_gui_json_output(gui_output_file)

            assert "command_line_opts" in json_data
            assert "config_opts" in json_data

            command_line_opts = json_data['command_line_opts']
            config_opts = json_data['config_opts']

            assert command_line_opts["input"] == runner.bam_file_name
            assert command_line_opts["bedfile"] == runner.bed_file_name

            assert config_opts["outputs"]["profiles"] is True
            assert config_opts["outputs"]["regions"] is True
            assert config_opts["pass"]["MINQCOV_MIN"] == 15
            assert config_opts["count_duplicate_reads"] is True
            assert config_opts["direction"] is True


if __name__ == "__main__":
    unittest.main()

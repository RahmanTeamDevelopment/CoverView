import os
import testutils.runners
import unittest


class TestCoverViewWithGuiOutput(unittest.TestCase):

    def test_coverview_runs_with_user_specified_gui_output_directory(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "gui": True,
                    "regions": True,
                    "profiles": True,
                    "gui_output_directory": os.path.join(os.getcwd(), "guiout")
                },

                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

    def test_coverview_runs_with_default_gui_output_file_names(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "gui": True,
                }
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0


if __name__ == "__main__":
    unittest.main()
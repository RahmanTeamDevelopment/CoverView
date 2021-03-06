import testutils.runners
import unittest


class TestCoverViewWithDirectionalOutput(unittest.TestCase):

    def test_coverview_runs_when_directional_output_is_requested(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            runner.add_config_data({
                "outputs": {
                    "profiles": True,
                },

                "count_duplicate_reads": True,
                "direction": True
            })

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0


if __name__ == "__main__":
    unittest.main()

import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverViewWithoutBedFile(unittest.TestCase):

    def test_produces_summary_output_when_run_with_no_bed_file(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            summary_output = testutils.output_checkers.load_coverview_summary_output(
                "output_summary.txt"
            )

            assert "Total" in summary_output
            assert "Unmapped" in summary_output
            assert "Mapped" in summary_output
            assert "1" in summary_output

            assert summary_output['Total']['RC'] == 1
            assert summary_output['Mapped']['RC'] == 1
            assert summary_output['Unmapped']['RC'] == 0

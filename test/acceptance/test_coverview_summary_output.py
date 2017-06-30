import testutils.runners
import testutils.output_checkers
import unittest


class TestCoverViewSummaryOutput(unittest.TestCase):

    def test_summary_output_with_zero_coverage_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            summary_output = testutils.output_checkers.load_coverview_summary_output(
                "output_summary.txt"
            )

            assert "Total" in summary_output
            assert "Unmapped" in summary_output
            assert "Mapped" in summary_output
            assert "1" in summary_output

            assert summary_output['Total']['RC'] == "0"
            assert summary_output['Total']['RCIN'] == "-"
            assert summary_output['Total']['RCOUT'] == "-"
            assert summary_output['Mapped']['RC'] == "0"
            assert summary_output['Mapped']['RCIN'] == "0"
            assert summary_output['Mapped']['RCOUT'] == "0"
            assert summary_output['Unmapped']['RC'] == "0"
            assert summary_output['Unmapped']['RCIN'] == "-"
            assert summary_output['Unmapped']['RCOUT'] == "-"

            chr_1 = summary_output['1']
            assert int(chr_1['RC']) == int(chr_1['RCIN']) + int(chr_1['RCOUT'])

    def test_summary_output_with_one_on_target_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            summary_output = testutils.output_checkers.load_coverview_summary_output(
                "output_summary.txt"
            )

            assert "Total" in summary_output
            assert "Unmapped" in summary_output
            assert "Mapped" in summary_output
            assert "1" in summary_output

            assert summary_output['Total']['RC'] == "1"
            assert summary_output['Total']['RCIN'] == "-"
            assert summary_output['Total']['RCOUT'] == "-"
            assert summary_output['Mapped']['RC'] == "1"
            assert summary_output['Mapped']['RCIN'] == "1"
            assert summary_output['Mapped']['RCOUT'] == "0"
            assert summary_output['Unmapped']['RC'] == "0"
            assert summary_output['Unmapped']['RCIN'] == "-"
            assert summary_output['Unmapped']['RCOUT'] == "-"

            chr_1 = summary_output['1']
            assert int(chr_1['RC']) == int(chr_1['RCIN']) + int(chr_1['RCOUT'])

    def test_summary_output_with_one_off_target_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 1000, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            summary_output = testutils.output_checkers.load_coverview_summary_output(
                "output_summary.txt"
            )

            assert "Total" in summary_output
            assert "Unmapped" in summary_output
            assert "Mapped" in summary_output
            assert "1" in summary_output

            assert summary_output['Total']['RC'] == "1"
            assert summary_output['Total']['RCIN'] == "-"
            assert summary_output['Total']['RCOUT'] == "-"
            assert summary_output['Mapped']['RC'] == "1"
            assert summary_output['Mapped']['RCIN'] == "0"
            assert summary_output['Mapped']['RCOUT'] == "1"
            assert summary_output['Unmapped']['RC'] == "0"
            assert summary_output['Unmapped']['RCIN'] == "-"
            assert summary_output['Unmapped']['RCOUT'] == "-"

            chr_1 = summary_output['1']
            assert int(chr_1['RC']) == int(chr_1['RCIN']) + int(chr_1['RCOUT'])

    def test_summary_output_with_one_on_target_and_one_off_target_read_in_bam(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 1))
            runner.add_reads(("1", 1000, 100, 1))
            runner.add_region(("1", 32, 132, "Region_1"))
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            summary_output = testutils.output_checkers.load_coverview_summary_output(
                "output_summary.txt"
            )

            assert "Total" in summary_output
            assert "Unmapped" in summary_output
            assert "Mapped" in summary_output
            assert "1" in summary_output

            assert summary_output['Total']['RC'] == "2"
            assert summary_output['Total']['RCIN'] == "-"
            assert summary_output['Total']['RCOUT'] == "-"
            assert summary_output['Mapped']['RC'] == "2"
            assert summary_output['Mapped']['RCIN'] == "1"
            assert summary_output['Mapped']['RCOUT'] == "1"
            assert summary_output['Unmapped']['RC'] == "0"
            assert summary_output['Unmapped']['RCIN'] == "-"
            assert summary_output['Unmapped']['RCOUT'] == "-"

            chr_1 = summary_output['1']
            assert int(chr_1['RC']) == int(chr_1['RCIN']) + int(chr_1['RCOUT'])

    def test_summary_output_reports_correct_number_of_unmapped_reads(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            read_length = 100
            num_unmapped_reads = 100
            num_mapped_reads = 50
            runner.add_reads(("1", 1000, read_length, num_mapped_reads))
            runner.add_unmapped_reads(read_length, num_unmapped_reads)
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            summary_output = testutils.output_checkers.load_coverview_summary_output(
                "output_summary.txt"
            )

            assert summary_output['Unmapped']['RC'] == str(num_unmapped_reads)
            assert summary_output['Mapped']['RC'] == str(num_mapped_reads)


if __name__ == "__main__":
    unittest.main()
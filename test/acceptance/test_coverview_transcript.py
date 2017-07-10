import testutils.runners
import testutils.output_checkers
import coverview.transcript
import unittest


class TestCoverViewWithTranscriptDatabase(unittest.TestCase):

    def test_runs_through_without_errors_when_transcript_db_is_used(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_transcript(
                coverview.transcript.Transcript(
                    "TEST_TRANSCRIPT_1",
                    "TEST_GENE_1",
                    "TEST_GENE_1",
                    "1",
                    "1",
                    50,
                    70,
                    55,
                    55,
                    70,
                    [
                        coverview.transcript.Exon(0, 55, 61),
                        coverview.transcript.Exon(1, 62, 70)
                    ]
                )
            )

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0


if __name__ == "__main__":
    unittest.main()

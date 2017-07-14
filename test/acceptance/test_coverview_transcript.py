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
                    ensembl_id="TEST_TRANSCRIPT_1",
                    gene_symbol="TEST_GENE_1",
                    gene_id="TEST_GENE_1",
                    chrom="1",
                    strand=1,
                    transcript_start=50,
                    transcript_end=70,
                    coding_start=0,
                    coding_start_genomic=50,
                    coding_end_genomic=65,
                    exons=[
                        coverview.transcript.Exon(0, 50, 60),
                        coverview.transcript.Exon(1, 62, 70)
                    ]
                )
            )
            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

    def test_transcript_coordinates_are_output_in_regions_file_when_transcript_regions_option_used(self):
        with testutils.runners.CoverViewTestRunner() as runner:
            runner.add_reads(("1", 32, 100, 0))
            runner.add_region(("1", 35, 40, "Region_1"))
            runner.add_region(("1", 50, 60, "Region_2"))
            runner.add_region(("1", 65, 70, "Region_3"))
            runner.add_config_data({
                "transcript": {
                    "regions": True
                }
            })

            runner.add_transcript(
                coverview.transcript.Transcript(
                    ensembl_id="TEST_TRANSCRIPT_1",
                    gene_symbol="TEST_GENE_1",
                    gene_id="TEST_GENE_1",
                    chrom="1",
                    strand=1,
                    transcript_start=50,
                    transcript_end=70,
                    coding_start=0,
                    coding_start_genomic=50,
                    coding_end_genomic=65,
                    exons=[
                        coverview.transcript.Exon(0, 50, 60),
                        coverview.transcript.Exon(1, 62, 70)
                    ]
                )
            )

            status_code = runner.run_coverview_and_get_exit_code()
            assert status_code == 0

            regions_output = testutils.output_checkers.load_coveriew_regions_output(
                "output_regions.txt"
            )

            assert regions_output['Region_1']["Start_transcript"] == ""
            assert regions_output['Region_1']["End_transcript"] == ""
            assert regions_output['Region_2']["Start_transcript"] == "TEST_GENE_1:TEST_TRANSCRIPT_1:c.1"
            assert regions_output['Region_2']["End_transcript"] == "TEST_GENE_1:TEST_TRANSCRIPT_1:c.10+1"
            assert regions_output['Region_3']["Start_transcript"] == "TEST_GENE_1:TEST_TRANSCRIPT_1:c.+1"
            assert regions_output['Region_3']["End_transcript"] == ""


if __name__ == "__main__":
    unittest.main()

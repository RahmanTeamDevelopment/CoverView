import coverview.transcript
import unittest


class TestCSNCalculationOnForwardTranscript(unittest.TestCase):

    def setUp(self):
        self.test_transcript = coverview.transcript.Transcript(
            ensembl_id="TEST_TRANSCRIPT_1",
            gene_symbol="TEST_GENE_1",
            gene_id="TEST_GENE_1",
            chrom="1",
            strand=1,
            transcript_start=0,
            transcript_end=100,
            coding_start=10,
            coding_start_genomic=10,
            coding_end_genomic=90,
            exons=[
                coverview.transcript.Exon(0, 0, 40),
                coverview.transcript.Exon(1, 50, 100)
                ]
        )

    def test_first_base_in_transcript(self):
        result = coverview.transcript.get_csn_coordinates(
            0,
            self.test_transcript
        )

        assert result == "c.-10"

    def test_first_coding_base_in_first_exon(self):
        result = coverview.transcript.get_csn_coordinates(
            10,
            self.test_transcript
        )

        assert result == "c.1"

    def test_last_base_in_first_exon(self):
        result = coverview.transcript.get_csn_coordinates(
            39,
            self.test_transcript
        )

        assert result == "c.30"

    def test_first_base_in_first_intron(self):
        result = coverview.transcript.get_csn_coordinates(
            40,
            self.test_transcript
        )

        assert result == "c.30+1"

    def test_middle_base_in_first_intron(self):
        result = coverview.transcript.get_csn_coordinates(
            44,
            self.test_transcript
        )

        assert result == "c.30+5"

    def test_last_base_in_first_intron(self):
        result = coverview.transcript.get_csn_coordinates(
            49,
            self.test_transcript
        )

        assert result == "c.31-1"

    def test_first_base_in_second_exon(self):
        result = coverview.transcript.get_csn_coordinates(
            50,
            self.test_transcript
        )

        assert result == "c.31"

    def test_last_base_in_second_exon(self):
        result = coverview.transcript.get_csn_coordinates(
            89,
            self.test_transcript
        )

        assert result == "c.70"

    def test_last_base_in_transcript(self):
        result = coverview.transcript.get_csn_coordinates(
            99,
            self.test_transcript
        )

        assert result == "c.+10"


if __name__ == "__main__":
    unittest.main()

import coverview.transcript
import unittest


class TestExon(unittest.TestCase):

    def test_exon_contains(self):
        start = 0
        end = 10
        exon = coverview.transcript.Exon(0, start, end)

        assert start - 1 not in exon
        assert start in exon
        assert end - 1 in exon
        assert end not in exon
        assert end + 1 not in exon

    def test_exon_distance_from(self):
        start = 0
        end = 10
        exon = coverview.transcript.Exon(0, start, end)

        assert exon.distance_from(start - 1) == -1
        assert exon.distance_from(start) == 0
        assert exon.distance_from(end -1) == 0
        assert exon.distance_from(end) == 1
        assert exon.distance_from(end + 1) == 2

    def test_get_coordinate_of_position_forward(self):
        start = 0
        end = 10
        exon = coverview.transcript.Exon(0, start, end)

        assert exon.get_coordinate_of_position_forward(start - 1) is None
        assert exon.get_coordinate_of_position_forward(start) == 0
        assert exon.get_coordinate_of_position_forward(start + 1) == 1
        assert exon.get_coordinate_of_position_forward(end - 1) == 9
        assert exon.get_coordinate_of_position_forward(end) is None
        assert exon.get_coordinate_of_position_forward(end + 1) is None

    def test_get_coordinate_of_position_reverse(self):
        start = 0
        end = 10
        exon = coverview.transcript.Exon(0, start, end)

        assert exon.get_coordinate_of_position_reverse(start - 1) is None
        assert exon.get_coordinate_of_position_reverse(start) == 9
        assert exon.get_coordinate_of_position_reverse(start + 1) == 8
        assert exon.get_coordinate_of_position_reverse(end - 1) == 0
        assert exon.get_coordinate_of_position_reverse(end) is None
        assert exon.get_coordinate_of_position_reverse(end + 1) is None


class TestTranscript(unittest.TestCase):
    def setUp(self):
        self.forward_transcript = coverview.transcript.Transcript(
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

        self.reverse_transcript = coverview.transcript.Transcript(
            ensembl_id="TEST_TRANSCRIPT_1",
            gene_symbol="TEST_GENE_1",
            gene_id="TEST_GENE_1",
            chrom="1",
            strand=-1,
            transcript_start=0,
            transcript_end=100,
            coding_start=10,
            coding_start_genomic=91,
            coding_end_genomic=5,
            exons=[
                coverview.transcript.Exon(1, 50, 100),
                coverview.transcript.Exon(0, 0, 40)
            ]
        )

    def test_is_position_in_utr_for_forward_transcript(self):
        assert self.forward_transcript.is_position_in_utr(0) is True
        assert self.forward_transcript.is_position_in_utr(9) is True
        assert self.forward_transcript.is_position_in_utr(10) is False
        assert self.forward_transcript.is_position_in_utr(89) is False
        assert self.forward_transcript.is_position_in_utr(90) is True
        assert self.forward_transcript.is_position_in_utr(100) is True
    
    def test_is_position_in_utr_for_reverse_transcript(self):
        assert self.reverse_transcript.is_position_in_utr(0) is True
        assert self.reverse_transcript.is_position_in_utr(5) is True
        assert self.reverse_transcript.is_position_in_utr(6) is False
        assert self.reverse_transcript.is_position_in_utr(91) is False
        assert self.reverse_transcript.is_position_in_utr(92) is True
        assert self.reverse_transcript.is_position_in_utr(100) is True

    def test_get_distance_from_coding_region(self):
        pass


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

    def _get_csn(self, position):
        return coverview.transcript.get_csn_coordinates(
            position,
            self.test_transcript
        )

    def test_first_base_in_transcript(self):
        assert self._get_csn(0) == "c.-10"

    def test_first_coding_base_in_first_exon(self):
        assert self._get_csn(10) == "c.1"

    def test_last_base_in_first_exon(self):
        assert self._get_csn(39) == "c.30"

    def test_first_base_in_first_intron(self):
        assert self._get_csn(40) == "c.30+1"

    def test_middle_base_in_first_intron(self):
        assert self._get_csn(44) == "c.30+5"

    def test_last_base_in_first_intron(self):
        assert self._get_csn(49) == "c.31-1"

    def test_first_base_in_second_exon(self):
        assert self._get_csn(50) == "c.31"

    def test_last_base_in_second_exon(self):
        assert self._get_csn(89) == "c.70"

    def test_last_base_in_transcript(self):
        assert self._get_csn(99) == "c.+10"


class TestCSNCalculationOnReverseTranscript(unittest.TestCase):

    def setUp(self):
        self.test_transcript = coverview.transcript.Transcript(
            ensembl_id="TEST_TRANSCRIPT_1",
            gene_symbol="TEST_GENE_1",
            gene_id="TEST_GENE_1",
            chrom="1",
            strand=-1,
            transcript_start=0,
            transcript_end=100,
            coding_start=10,
            coding_start_genomic=91,
            coding_end_genomic=5,
            exons=[
                coverview.transcript.Exon(1, 50, 100),
                coverview.transcript.Exon(0, 0, 40)
                ]
        )

    def _get_csn(self, position):
        return coverview.transcript.get_csn_coordinates(
            position,
            self.test_transcript
        )

    def test_first_base_in_transcript(self):
        assert self._get_csn(99) == "c.-8"

    def test_first_coding_base_in_first_exon(self):
        assert self._get_csn(91) == "c.1"
    #
    # def test_last_base_in_first_exon(self):
    #     assert self._get_csn(50) == "c.50"
    #
    # def test_first_base_in_first_intron(self):
    #     assert self._get_csn(49) == "c.50+1"
    #
    # def test_middle_base_in_first_intron(self):
    #     assert self._get_csn(45) == "c.50+5"
    #
    # def test_last_base_in_first_intron(self):
    #     assert self._get_csn(40) == "c.51-1"
    #
    # def test_first_base_in_second_exon(self):
    #     assert self._get_csn(39) == "c.51"
    #
    # def test_last_coding_base_in_second_exon(self):
    #     assert self._get_csn(6) == "c.84"
    #
    # def test_last_base_in_transcript(self):
    #     assert self._get_csn(0) == "c.+5"

if __name__ == "__main__":
    unittest.main()

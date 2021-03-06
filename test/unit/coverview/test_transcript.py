import coverview_.transcript
import unittest


class TestExon(unittest.TestCase):

    def test_exon_contains(self):
        start = 0
        end = 10
        exon = coverview_.transcript.Exon(0, start, end)

        assert start - 1 not in exon
        assert start in exon
        assert end - 1 in exon
        assert end not in exon
        assert end + 1 not in exon

    def test_exon_distance_from(self):
        start = 0
        end = 10
        exon = coverview_.transcript.Exon(0, start, end)

        assert exon.distance_from(start - 1) == -1
        assert exon.distance_from(start) == 0
        assert exon.distance_from(end - 1) == 0
        assert exon.distance_from(end) == 1
        assert exon.distance_from(end + 1) == 2

    def test_get_coordinate_of_position_forward(self):
        start = 0
        end = 10
        exon = coverview_.transcript.Exon(0, start, end)

        assert exon.get_coordinate_of_position_forward(start - 1) is None
        assert exon.get_coordinate_of_position_forward(start) == 0
        assert exon.get_coordinate_of_position_forward(start + 1) == 1
        assert exon.get_coordinate_of_position_forward(end - 1) == 9
        assert exon.get_coordinate_of_position_forward(end) is None
        assert exon.get_coordinate_of_position_forward(end + 1) is None

    def test_get_coordinate_of_position_reverse(self):
        start = 0
        end = 10
        exon = coverview_.transcript.Exon(0, start, end)

        assert exon.get_coordinate_of_position_reverse(start - 1) is None
        assert exon.get_coordinate_of_position_reverse(start) == 9
        assert exon.get_coordinate_of_position_reverse(start + 1) == 8
        assert exon.get_coordinate_of_position_reverse(end - 1) == 0
        assert exon.get_coordinate_of_position_reverse(end) is None
        assert exon.get_coordinate_of_position_reverse(end + 1) is None


class TestTranscript(unittest.TestCase):
    def setUp(self):
        self.forward_transcript = coverview_.transcript.Transcript(
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
                coverview_.transcript.Exon(0, 0, 40),
                coverview_.transcript.Exon(1, 50, 100)
            ]
        )

        self.reverse_transcript = coverview_.transcript.Transcript(
            ensembl_id="TEST_TRANSCRIPT_1",
            gene_symbol="TEST_GENE_1",
            gene_id="TEST_GENE_1",
            chrom="1",
            strand=-1,
            transcript_start=0,
            transcript_end=100,
            coding_start=10,
            coding_start_genomic=89,
            coding_end_genomic=5,
            exons=[
                coverview_.transcript.Exon(1, 50, 100),
                coverview_.transcript.Exon(0, 0, 40)
            ]
        )

    def test_total_length_of_coding_sequence(self):
        assert self.forward_transcript.total_length_of_coding_sequence == 70
        assert self.reverse_transcript.total_length_of_coding_sequence == 74


class TestCreateTranscriptFromDataBaseTextLine(unittest.TestCase):

    def test_create_forward_transcript_from_database_line(self):
        line = "ENST00000379389\tISG15\tENSG00000187608\t+/1.1kb/2/0.7kb/165" \
               "\t1\t1\t948802\t949920\t152\t948954\t949858\t948802\t948956\t949363\t949920"

        transcript = coverview_.transcript.create_transcript_from_line_of_old_database(
            line
        )

        assert transcript.ensembl_id == "ENST00000379389"
        assert transcript.gene_symbol == "ISG15"
        assert transcript.gene_id == "ENSG00000187608"
        assert transcript.chrom == "1"
        assert transcript.strand == 1
        assert transcript.transcript_start == 948802
        assert transcript.transcript_end == 949920
        assert transcript.coding_start == 151
        assert transcript.coding_start_genomic == 948953
        assert transcript.coding_start_genomic == transcript.transcript_start + transcript.coding_start
        assert transcript.coding_end_genomic == 949858
        assert len(transcript.exons) == 2
        assert transcript.exons[0].start == 948802
        assert transcript.exons[0].end == 948956
        assert transcript.exons[1].start == 949363
        assert transcript.exons[1].end == 949920

    def test_create_reverse_transcript_from_database_line(self):
        line = "ENST00000422725\tC1orf233\tENSG00000228594\t-/2.1kb/1/2.1kb/226" \
               "\t1\t-1\t1533391\t1535476\t82\t1535395\t1534715\t1533391\t1535476"

        transcript = coverview_.transcript.create_transcript_from_line_of_old_database(
            line
        )

        assert transcript.ensembl_id == "ENST00000422725"
        assert transcript.gene_symbol == "C1orf233"
        assert transcript.gene_id == "ENSG00000228594"
        assert transcript.chrom == "1"
        assert transcript.strand == -1
        assert transcript.transcript_start == 1533391
        assert transcript.transcript_end == 1535476
        assert transcript.coding_start == 81
        assert transcript.coding_start_genomic == 1535394
        assert transcript.coding_start_genomic == transcript.transcript_end - 1 - transcript.coding_start
        assert transcript.coding_end_genomic == 1534713
        assert len(transcript.exons) == 1
        assert transcript.exons[0].start == 1533391
        assert transcript.exons[0].end == 1535476


class TestCSNCalculationOnForwardTranscript(unittest.TestCase):

    def setUp(self):
        self.test_transcript = coverview_.transcript.Transcript(
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
                coverview_.transcript.Exon(0, 0, 40),
                coverview_.transcript.Exon(1, 50, 100)
                ]
        )

    def _get_csn(self, position):
        return coverview_.transcript.get_csn_coordinates(
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
        self.test_transcript = coverview_.transcript.Transcript(
            ensembl_id="TEST_TRANSCRIPT_1",
            gene_symbol="TEST_GENE_1",
            gene_id="TEST_GENE_1",
            chrom="1",
            strand=-1,
            transcript_start=0,
            transcript_end=100,
            coding_start=10,
            coding_start_genomic=89,
            coding_end_genomic=5,
            exons=[
                coverview_.transcript.Exon(1, 50, 100),
                coverview_.transcript.Exon(0, 0, 40)
                ]
        )

    def _get_csn(self, position):
        return coverview_.transcript.get_csn_coordinates(
            position,
            self.test_transcript
        )

    def test_first_base_in_transcript(self):
        assert self._get_csn(99) == "c.-10"

    def test_first_coding_base_in_first_exon(self):
        assert self._get_csn(89) == "c.1"

    def test_last_base_in_first_exon(self):
        assert self._get_csn(50) == "c.40"

    def test_first_base_in_first_intron(self):
        assert self._get_csn(49) == "c.40+1"

    def test_middle_base_in_first_intron(self):
        assert self._get_csn(45) == "c.40+5"

    def test_last_base_in_first_intron(self):
        assert self._get_csn(40) == "c.41-1"

    def test_first_base_in_second_exon(self):
        assert self._get_csn(39) == "c.41"

    def test_last_coding_base_in_second_exon(self):
        assert self._get_csn(6) == "c.74"

    def test_last_base_in_transcript(self):
        assert self._get_csn(0) == "c.+6"


class TestCSNCalculationOnForwardTranscriptWhenCodingStartIsInSecondExon(unittest.TestCase):

    def setUp(self):
        self.test_transcript = coverview_.transcript.Transcript(
            ensembl_id="TEST_TRANSCRIPT_1",
            gene_symbol="TEST_GENE_1",
            gene_id="TEST_GENE_1",
            chrom="1",
            strand=1,
            transcript_start=0,
            transcript_end=100,
            coding_start=60,
            coding_start_genomic=70,
            coding_end_genomic=90,
            exons=[
                coverview_.transcript.Exon(0, 0, 40),
                coverview_.transcript.Exon(1, 50, 100)
                ]
        )

    def _get_csn(self, position):
        return coverview_.transcript.get_csn_coordinates(
            position,
            self.test_transcript
        )

    def test_first_base_in_transcript(self):
        assert self._get_csn(0) == "c.-60"

    def test_first_coding_base(self):
        assert self._get_csn(70) == "c.1"

    def test_last_base_in_first_exon(self):
        assert self._get_csn(39) == "c.-21"

    def test_first_base_in_first_intron(self):
        assert self._get_csn(40) == "c.-21+1"

    def test_last_base_in_transcript(self):
        assert self._get_csn(99) == "c.+10"


class TestCSNForBRCA2(unittest.TestCase):
    def setUp(self):
        line = "ENST00000380152\tBRCA2\tENSG00000139618\t+/83.7kb/27/10.9kb/3418\t13\t1" \
               "\t32889610\t32973347\t234\t32890598\t32972907\t32889610\t32889804\t32890558\t32890664" \
               "\t32893213\t32893462\t32899212\t32899321\t32900237\t32900287\t32900378\t32900419\t32900635" \
               "\t32900750\t32903579\t32903629\t32905055\t32905167\t32906408\t32907524\t32910401\t32915333" \
               "\t32918694\t32918790\t32920963\t32921033\t32928997\t32929425\t32930564\t32930746\t32931878" \
               "\t32932066\t32936659\t32936830\t32937315\t32937670\t32944538\t32944694\t32945092\t32945237" \
               "\t32950806\t32950928\t32953453\t32953652\t32953886\t32954050\t32954143\t32954282\t32968825" \
               "\t32969070\t32971034\t32971181\t32972298\t32973347"

        self.brca2 = coverview_.transcript.create_transcript_from_line_of_old_database(
            line
        )

    def _get_csn(self, position):
        return coverview_.transcript.get_csn_coordinates(
            position,
            self.brca2
        )

    def test_csn_for_brca2(self):
        assert self.brca2.coding_start == 233
        assert self.brca2.coding_start_genomic == 32890597
        assert self.brca2.coding_end_genomic == 32972907
        assert self.brca2.strand == 1
        assert len(self.brca2.exons) == 27

        for exon in self.brca2.exons:
            assert exon.length > 0

        assert self._get_csn(32889610) == "c.-233"
        assert self._get_csn(32890587) == "c.-10"
        assert self._get_csn(32890588) == "c.-9"
        assert self._get_csn(32890589) == "c.-8"
        assert self._get_csn(32890590) == "c.-7"
        assert self._get_csn(32890591) == "c.-6"
        assert self._get_csn(32890592) == "c.-5"
        assert self._get_csn(32890593) == "c.-4"
        assert self._get_csn(32890594) == "c.-3"
        assert self._get_csn(32890595) == "c.-2"
        assert self._get_csn(32890596) == "c.-1"
        assert self._get_csn(32890597) == "c.1"
        assert self._get_csn(32890598) == "c.2"
        assert self._get_csn(32972907) == "c.+1"
        assert self._get_csn(32972908) == "c.+2"
        assert self._get_csn(32972909) == "c.+3"


if __name__ == "__main__":
    unittest.main()

import bamgen
import os
import pysam
import unittest
import uuid


class TestMockFastaFile(unittest.TestCase):

    def test_returns_A_for_single_base_fetch(self):
        fasta_file = bamgen.MockReferenceFile()
        assert fasta_file.fetch("1", 0, 1) == "A"

    def test_returns_AA_for_two_base_fetch(self):
        fasta_file = bamgen.MockReferenceFile()
        assert fasta_file.fetch("1", 0, 2) == "AA"


class TestSimpleBamFileGeneration(unittest.TestCase):
    """
    Test generating BAM files, and reading data from those BAM files. Each
    BAM file has a unique file-name, which is a hash generated through the
    uuid module. The BAM files and their indexes are deleted after each test.
    """
    def setUp(self):
        self.unique_bam_file_name = str(uuid.uuid4())
        self.unique_index_file_name = self.unique_bam_file_name + ".bai"

    def tearDown(self):
        os.remove(self.unique_bam_file_name)
        os.remove(self.unique_index_file_name)

    def test_can_generate_empty_bam_file(self):
        ref_file = bamgen.MockReferenceFile()
        regions = [
            ("1", 32, 100, 0)
        ]

        bamgen.generate_bam_file(
            self.unique_bam_file_name,
            ref_file,
            regions
        )

        with pysam.AlignmentFile(self.unique_bam_file_name, 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 0

    def test_can_generate_bam_file_with_single_read(self):
        ref_file = bamgen.MockReferenceFile()
        regions = [
            ("1", 32, 100, 1)
        ]

        bamgen.generate_bam_file(
            self.unique_bam_file_name,
            ref_file,
            regions
        )

        with pysam.AlignmentFile(self.unique_bam_file_name, 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 1

    def test_can_generate_bam_file_with_two_reads(self):
        ref_file = bamgen.MockReferenceFile()
        regions = [
            ("1", 32, 100, 2)
        ]

        bamgen.generate_bam_file(
            self.unique_bam_file_name,
            ref_file,
            regions
        )

        with pysam.AlignmentFile(self.unique_bam_file_name, 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 2

    def test_can_generate_bam_file_with_many_reads(self):
        ref_file = bamgen.MockReferenceFile()
        regions = [
            ("1", 32, 100, 1000)
        ]

        bamgen.generate_bam_file(
            self.unique_bam_file_name,
            ref_file,
            regions
        )

        with pysam.AlignmentFile(self.unique_bam_file_name, 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 1000


if __name__ == "__main__":
    unittest.main()
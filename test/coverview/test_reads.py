import bamgen
import coverview
import os
import pysam
import unittest
import uuid


class TestReadArray(unittest.TestCase):
    """
    Here we are testing the ability to read a large set of reads into
    memory using the ReadArray class, and to query sub-regions of the
    data.
    """
    def setUp(self):
        self.unique_bam_file_name = str(uuid.uuid4())
        self.unique_index_file_name = self.unique_bam_file_name + ".bai"

    def tearDown(self):
        os.remove(self.unique_bam_file_name)
        os.remove(self.unique_index_file_name)

    def test_empty_read_array_counts_zero_reads_in_interval(self):
        references = ["1"]
        reference_lengths = [1000]
        ref_file = bamgen.MockReferenceFile(references, reference_lengths)
        regions = [
            ("1", 32, 100, 0)
        ]

        bamgen.generate_bam_file(
            self.unique_bam_file_name,
            ref_file,
            regions
        )

        read_array = coverview.reads.pyReadArray()

        with pysam.AlignmentFile(self.unique_bam_file_name, 'rb') as bam_file:
            for read in bam_file:
                read_array.append(read)

        assert read_array.count_reads_in_interval(0, 1000) == 0

    def test_read_array_with_single_counts_one_read_for_interval_overlapping_read(self):
        references = ["1"]
        reference_lengths = [1000]
        ref_file = bamgen.MockReferenceFile(references, reference_lengths)
        regions = [
            ("1", 32, 100, 1)
        ]

        bamgen.generate_bam_file(
            self.unique_bam_file_name,
            ref_file,
            regions
        )

        read_array = coverview.reads.pyReadArray()

        with pysam.AlignmentFile(self.unique_bam_file_name, 'rb') as bam_file:
            for read in bam_file:
                read_array.append(read)

        assert read_array.count_reads_in_interval(10, 50) == 1


if __name__ == "__main__":
    unittest.main()
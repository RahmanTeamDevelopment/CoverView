import bamgen.bamgen
import coverview.statistics
import coverview.reads
import os
import pysam
import unittest
import uuid


def load_bam_into_read_array(file_name):
    """
    Utility function for creating a read array from the contents of a sorted BAM file    
    """
    read_array = coverview.reads.pyReadArray()

    with pysam.AlignmentFile(file_name, 'rb') as bam_file:
        for read in bam_file:
            read_array.append(read)

    return read_array


class TestReadArray(unittest.TestCase):
    """
    Here we are testing the ability to read a set of reads from a BAM file into
    memory using the ReadArray class, and to query sub-regions of the
    reads.
    """
    def setUp(self):
        self.unique_bam_file_name = str(uuid.uuid4())
        self.unique_index_file_name = self.unique_bam_file_name + ".bai"

    def tearDown(self):
        os.remove(self.unique_bam_file_name)
        os.remove(self.unique_index_file_name)

    def test_empty_read_array_counts_zero_reads_in_interval(self):
        read_sets = [
            ("1", 32, 100, 0)
        ]

        bamgen.bamgen.make_bam_file(self.unique_bam_file_name, read_sets)
        read_array = load_bam_into_read_array(self.unique_bam_file_name)

        assert read_array.count_reads_in_interval(32, 132) == 0

    def test_read_array_with_single_read_counts_one_for_interval_overlapping_read(self):
        read_sets = [
            ("1", 32, 100, 1)
        ]

        bamgen.bamgen.make_bam_file(self.unique_bam_file_name, read_sets)
        read_array = load_bam_into_read_array(self.unique_bam_file_name)

        assert read_array.count_reads_in_interval(32, 132) == 1

    def test_read_array_with_single_read_counts_zero_for_interval_not_overlapping_read(self):
        read_sets = [
            ("1", 32, 100, 1)
        ]

        bamgen.bamgen.make_bam_file(self.unique_bam_file_name, read_sets)
        read_array = load_bam_into_read_array(self.unique_bam_file_name)

        assert read_array.count_reads_in_interval(0, 31) == 0
        assert read_array.count_reads_in_interval(132, 133) == 0


if __name__ == "__main__":
    unittest.main()
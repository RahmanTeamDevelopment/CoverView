import bamgen
import coverview
import os
import pysam
import unittest
import uuid


def make_bam_file(file_name, read_sets):
    """
    Utility function to create a BAM file with the specified configuration. We create a mock
    reference with the union of all references in the read_sets and each reference sequence is
    at least as long as the max(2 * (start_position + read_length)), which is longer than we need but
    that's ok.
    """
    references = sorted(
        set([x[0] for x in read_sets])
    )

    num_read_sets = len(read_sets)
    longest_ref = max(2* (x[1] + x[2]) for x in read_sets)
    reference_lengths = [longest_ref] * num_read_sets
    ref_file = bamgen.MockReferenceFile(references, reference_lengths)

    bamgen.generate_bam_file(
        file_name,
        ref_file,
        read_sets
    )


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

        make_bam_file(self.unique_bam_file_name, read_sets)
        read_array = load_bam_into_read_array(self.unique_bam_file_name)

        assert read_array.count_reads_in_interval(32, 132) == 0

    def test_read_array_with_single_read_counts_one_for_interval_overlapping_read(self):
        read_sets = [
            ("1", 32, 100, 1)
        ]

        make_bam_file(self.unique_bam_file_name, read_sets)
        read_array = load_bam_into_read_array(self.unique_bam_file_name)

        assert read_array.count_reads_in_interval(32, 132) == 1

    def test_read_array_with_single_read_counts_zero_for_interval_not_overlapping_read(self):
        read_sets = [
            ("1", 32, 100, 1)
        ]

        make_bam_file(self.unique_bam_file_name, read_sets)
        read_array = load_bam_into_read_array(self.unique_bam_file_name)

        assert read_array.count_reads_in_interval(0, 31) == 0
        assert read_array.count_reads_in_interval(132, 133) == 0


if __name__ == "__main__":
    unittest.main()
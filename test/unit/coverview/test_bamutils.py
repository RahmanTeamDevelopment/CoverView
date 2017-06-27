import coverview.bamutils
import unittest


class TestBamIndexStats(unittest.TestCase):

    def test_load_bam_index_stats_from_text_representation(self):
        chromosome = "chr1"
        chromosome_length = 1000000
        num_mapped_reads = 100
        num_unmapped_reads = 1000

        index_stats_text = "{}\t{}\t{}\t{}".format(
            chromosome,
            chromosome_length,
            num_mapped_reads,
            num_unmapped_reads
        )

        index_stats = coverview.bamutils.extract_bam_index_stats(
            [index_stats_text]
        )

        assert chromosome in index_stats.data_by_chrom
        assert index_stats.get_length_of_chromosome(chromosome) == chromosome_length
        assert index_stats.get_num_mapped_reads_for_chromosome(chromosome) == num_mapped_reads
        assert index_stats.get_num_unmapped_reads_for_chromosome(chromosome) == num_unmapped_reads
        assert index_stats.get_num_total_reads_for_chromosome(chromosome) == num_unmapped_reads + num_mapped_reads

    def test_raises_standard_error_if_invalid_chromosome_name_passed(self):
        chromosome = "chr1"
        chromosome_length = 1000000
        num_mapped_reads = 100
        num_unmapped_reads = 1000

        index_stats_text = "{}\t1{}\t{}\t{}\n".format(
            chromosome,
            chromosome_length,
            num_mapped_reads,
            num_unmapped_reads
        )

        index_stats = coverview.bamutils.extract_bam_index_stats(
            [index_stats_text]
        )

        self.assertRaises(
            StandardError,
            index_stats.get_num_mapped_reads_for_chromosome,
            "chr2"
        )

    def test_counting_total_reads_in_bam_file(self):
        index_stats = coverview.bamutils.extract_bam_index_stats([
            "chr1\t1000000\t1000\t2000",
            "chr2\t1000000\t5000\t3000",
            "*\t0\t0\t120000"
        ])

        assert index_stats.get_total_reads_in_bam() == 131000

    def test_counting_total_mapped_reads_in_bam_file(self):
        index_stats = coverview.bamutils.extract_bam_index_stats([
            "chr1\t1000000\t1000\t2000",
            "chr2\t1000000\t5000\t3000",
            "*\t0\t0\t120000"
        ])

        assert index_stats.get_total_mapped_reads_in_bam() == 6000

    def test_counting_total_unmapped_reads_in_bam_file(self):
        index_stats = coverview.bamutils.extract_bam_index_stats([
            "chr1\t1000000\t1000\t2000",
            "chr2\t1000000\t5000\t3000",
            "*\t0\t0\t120000"
        ])

        assert index_stats.get_total_unmapped_reads_in_bam() == 125000


if __name__ == "__main__":
    unittest.main()
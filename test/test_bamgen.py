import bamgen
import os
import pysam
import unittest


class TestSimpleBamFileGeneration(unittest.TestCase):

    def tearDown(self):
        os.remove("test1.bam")
        os.remove("test1.bam.bai")
        
    def test_can_generate_bam_file_with_single_read(self):
        config = {
            "bam_files":
                [
                    {
                        "file_name": "test1.bam",
                        "reference": "/Users/arimmer/Work/Reference/Reference/human_g1k_v37.fa",
                        "reads":
                        [
                            {
                                "chrom": "1",
                                "start_pos": 32,
                                "read_length": 100,
                                "num_reads": 1
                            }
                        ]
                    }
                ]
        }

        bamgen.generate_bam_files(config)

        with pysam.AlignmentFile("test1.bam", 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 1

    def test_can_generate_bam_file_with_two_reads(self):
        config = {
            "bam_files":
                [
                    {
                        "file_name": "test1.bam",
                        "reference": "/Users/arimmer/Work/Reference/Reference/human_g1k_v37.fa",
                        "reads":
                        [
                            {
                                "chrom": "1",
                                "start_pos": 32,
                                "read_length": 100,
                                "num_reads": 2
                            }
                        ]
                    }
                ]
        }

        bamgen.generate_bam_files(config)

        with pysam.AlignmentFile("test1.bam", 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 2

    def test_can_generate_bam_file_with_many_reads(self):
        config = {
            "bam_files":
                [
                    {
                        "file_name": "test1.bam",
                        "reference": "/Users/arimmer/Work/Reference/Reference/human_g1k_v37.fa",
                        "reads":
                        [
                            {
                                "chrom": "1",
                                "start_pos": 32,
                                "read_length": 100,
                                "num_reads": 1000
                            }
                        ]
                    }
                ]
        }

        bamgen.generate_bam_files(config)

        with pysam.AlignmentFile("test1.bam", 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 1000




if __name__ == "__main__":
    unittest.main()
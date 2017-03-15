import bamgen
import pysam
import unittest


class TestBamGen(unittest.TestCase):
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

        with pysam.AlignmentFile(config['file_name'], 'rb') as bam_file:
            assert bam_file.count("1", 0, 100) == 1
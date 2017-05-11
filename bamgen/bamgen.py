import array
import json
import pysam
import sys


class MockReferenceFile(object):
    """
    This class is used in place of a real FASTA reference file, when we want to
    generate test BAM files, but don't care too much about the reads having a
    realistic sequence, and also don't want to have a large reference FASTA file
    lying around.
    """
    def __init__(self, references, lengths):
        self.references = references
        self.lengths = lengths

    def fetch(self, chrom, start_pos, end_pos):
        assert chrom in self.references
        return 'A' * (end_pos - start_pos)


class PerfectReadGenerator(object):
    """
    Instances of this class can be used to generate 'perfect' reads. These
    reads will all have sequences matching the reference exactly, and will have
    no alignment gaps. All bases will be high quality, mapping quality will be
    high, and the reads will not be marked as duplicates.
    """
    HIGH_QUAL = 60  # Used for all base and mapping qualities

    def __init__(self, ref_file):
        """
        """
        self.ref = ref_file

    def generate_unpaired_reads(self, chrom, chrom_id, pos, length, n,
                                read_id=0):
        """
        Generate n perfect (i.e. no mis-matches or alignment gaps, not
        duplicates, all bases are high quality, and the mapping quality
        is high) un-paired reads.
        """
        seq = self.ref.fetch(chrom, pos, pos+length)
        qual = array.array('b', [self.HIGH_QUAL]*len(seq))
        cigar = ((0, len(seq)),)  # Perfect match

        for i in xrange(n):
            read = pysam.AlignedRead()
            read.query_name = "simulated_read_{}".format(read_id + i)
            read.query_sequence = seq
            read.query_qualities = qual
            read.reference_id = chrom_id
            read.reference_start = pos
            read.mapping_quality = self.HIGH_QUAL
            read.cigar = cigar
            read.next_reference_id = -1
            read.next_reference_start = -1
            read.is_unmapped = 0
            read.is_duplicate = 0
            read.is_paired = 0
            yield read


def generate_bam_files(config):
    for bam_file_config in config['bam_files']:
        ref_file_name = bam_file_config['reference']
        ref_file = pysam.FastaFile(ref_file_name)
        bam_file_name = bam_file_config['file_name']
        read_gen = PerfectReadGenerator(ref_file)

        with pysam.AlignmentFile(
            bam_file_name,
            'wb',
            reference_names=ref_file.references,
            reference_lengths=ref_file.lengths
        ) as bam_file:

            regions = bam_file_config['reads']

            for region in regions:
                chrom = region['chrom']
                start_pos = region['start_pos']
                read_length = region['read_length']
                num_reads = region['num_reads']
                chrom_id = bam_file.get_tid(chrom)

                for read in read_gen.generate_unpaired_reads(
                        chrom,
                        chrom_id,
                        start_pos,
                        read_length,
                        num_reads,
                        read_id=0
                ):
                    bam_file.write(read)

        pysam.index(bam_file_name)


if __name__ == "__main__":
    config_file_name = sys.argv[1]

    with open(config_file_name, 'r') as config_file:
        config = json.load(config_file)
        generate_bam_files(config)

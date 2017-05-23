import array
import pysam
import uuid


_high_quality = 60


def make_bam_file(file_name, read_sets):
    """
    Utility function to create a BAM file with the specified configuration. We create a mock
    reference with the union of all references in the read_sets and each reference sequence is
    at least as long as the max(2 * (start_position + read_length)), which is longer than we need but
    that's ok.
    """
    ref_file = MockReferenceFile()

    generate_bam_file(
        file_name,
        ref_file,
        read_sets
    )


class MockReferenceFile(object):
    """
    This class is used in place of a real FASTA reference file, when we want to
    generate test BAM files, but don't care too much about the reads having a
    realistic sequence, and also don't want to have a large reference FASTA file
    lying around.
    """
    def __init__(self):
        pass

    def fetch(self, chrom, start_pos, end_pos):
        return 'A' * (end_pos - start_pos)

    def getReferenceLength(self, chrom):
        """
        The name of this function matches the one in pysam.Fastafile
        """
        return 2**64  # A very, very large number.


def create_unpaired_read(
        sequence,
        qualities,
        reference_id,
        start_position,
        mapping_quality,
        cigar,
        is_unmapped,
        is_duplicate
):
    read = pysam.AlignedRead()
    read.query_name = "simulated_read_{}".format(
        str(uuid.uuid4())
    )
    read.query_sequence = sequence
    read.query_qualities = qualities
    read.reference_id = reference_id
    read.reference_start = start_position
    read.mapping_quality = mapping_quality
    read.cigar = cigar
    read.next_reference_id = -1
    read.next_reference_start = -1
    read.is_unmapped = is_unmapped
    read.is_duplicate = is_duplicate
    read.is_paired = 0
    return read


def create_perfect_unpaired_read(
        reference_file,
        reference_id,
        chromosome,
        start_position,
        read_length
):
    mapping_quality = _high_quality
    sequence = reference_file.fetch(chromosome, start_position, start_position + read_length)
    qualities = array.array('b', [_high_quality] * read_length)
    cigar = ((0, read_length),)  # Perfect match

    return create_unpaired_read(
        sequence,
        qualities,
        reference_id,
        start_position,
        mapping_quality,
        cigar,
        is_unmapped=0,
        is_duplicate=0
    )


def generate_bam_file(bam_file_name, reference_file, regions):
    """
    Create a new BAM file and write reads to that file. The number and details
    of the reads are specified in the 'regions' parameter.
    """
    references = sorted(set(x[0] for x in regions))
    lengths = [max(2 * (x[1] + x[2]) for x in regions)] * len(references)

    with pysam.AlignmentFile(
        bam_file_name,
        'wb',
        reference_names=references,
        reference_lengths=lengths
    ) as bam_file:

        for chromosome, start_position, read_length, num_reads in regions:
            reference_id = bam_file.get_tid(chromosome)

            for read_index in xrange(num_reads):
                new_read = create_perfect_unpaired_read(
                    reference_file,
                    reference_id,
                    chromosome,
                    start_position,
                    read_length
                )

                bam_file.write(new_read)

    pysam.index(bam_file_name)

import array
import pysam
import sys


HIGH_QUAL = 100


def write_perfect_unpaired_read_to_bam(bam, chrom, pos, seq, read_id):
    """
    Write a single perfect (i.e. no mis-matches or alignment gaps, is not a
    duplicate, all bases are high quality, and the mapping quality is high)
    un-paired read to the specified BAM file.
    """
    read = pysam.AlignedRead()

    read.query_name = "simulated_read_{}".format(read_id)
    read.query_sequence = seq
    read.query_qualities = array.array('b', [HIGH_QUAL]*len(seq))
    read.reference_id = bam.get_tid(chrom)
    read.reference_start = pos
    read.mapping_quality = HIGH_QUAL
    read.cigar = ((0, len(seq)),)  # Perfect match
    read.next_reference_id = -1
    read.next_reference_start = -1
    read.is_unmapped = 0
    read.is_duplicate = 0
    read.is_paired = 0

    bam.write(read)


def write_n_perfect_unpaired_reads_to_bam(bam, ref, chrom, pos, length,
                                          read_id, n):
    """
    Write a number of identical, perfect (i.e. no mis-matches or alignment
    gaps, is not a duplicate, all bases are high quality, and the mapping
    quality is high) un-paired reads to the specified BAM file.
    """
    seq = ref.fetch(chrom, pos, pos+length)

    for i in xrange(n):
        write_perfect_unpaired_read_to_bam(bam, chrom, pos, seq, read_id + i)


if __name__ == "__main__":
    ref_file = pysam.FastaFile(sys.argv[1])

    bam_file = pysam.AlignmentFile(
        "test.bam",
        'wb',
        reference_names=ref_file.references,
        reference_lengths=ref_file.lengths
    )

    seq = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    chrom = "1"
    start_pos = 32
    read_length = 100

    write_n_perfect_unpaired_reads_to_bam(
        bam_file,
        ref_file,
        chrom,
        start_pos,
        read_length,
        0,
        10000)

    bam_file.close()

import pysam
import sys


def write_read_to_bam(bam, seq, qual, chrom, pos, map_qual, cigar, read_id):
    """
    Write a single read to the specified BAM file.
    """
    read = pysam.AlignedRead()

    read.query_name = "simulated_read_{}".format(read_id)
    read.query_sequence = seq
    read.query_qualities = pysam.qualitystring_to_array(qual)
    read.reference_id = chrom
    read.reference_start = pos
    read.mapping_quality = map_qual
    read.cigar = cigar
    read.next_reference_id = -1
    read.next_reference_start = -1
    read.is_unmapped = 0
    read.is_duplicate = 0
    read.is_paired = 0

    bam.write(read)


if __name__ == "__main__":
    ref_file = pysam.Fastafile(sys.argv[1])

    bam_file = pysam.Samfile(
        "test.bam",
        'wb',
        reference_names=ref_file.references,
        reference_lengths=ref_file.lengths
    )

    for i in xrange(10000):
        seq = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        qual = "<"*len(seq)
        tid = 0  # chr1
        start_pos = 32
        map_qual = 50
        cigar = ((0, 10),)

        write_read_to_bam(
            bam_file,
            seq,
            qual,
            tid,
            start_pos,
            map_qual,
            cigar,
            i)

    bam_file.close()

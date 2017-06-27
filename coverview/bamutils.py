import pysam


class BamIndexStatsRow(object):
    def __init__(
            self,
            chromosome,
            chromosome_length,
            num_reads,
            num_mapped_reads,
            num_unmapped_reads
    ):
        self.chromosome = chromosome
        self.chromosome_length = chromosome_length
        self.num_reads = num_reads
        self.num_mapped_reads = num_mapped_reads
        self.num_unmapped_reads = num_unmapped_reads


class BamIndexStats(object):
    """
    Utility class for wrapping the data returned from pysam.idxstats
    """
    def __init__(self):
        self.rows = []
        self.data_by_chrom = {}

    def add_row(self, row):
        self.rows.append(row)
        self.data_by_chrom[row.chromosome] = row

    def _get_stats_for_chromosome(self, chromosome):
        if chromosome not in self.data_by_chrom:
            raise StandardError(
                "Invalid chromosome name ({})".format(chromosome)
            )
        else:
            return self.data_by_chrom[chromosome]

    def get_num_mapped_reads_for_chromosome(self, chromosome):
        stats = self._get_stats_for_chromosome(chromosome)
        return stats.num_mapped_reads

    def get_num_unmapped_reads_for_chromosome(self, chromosome):
        stats = self._get_stats_for_chromosome(chromosome)
        return stats.num_mapped_reads

    def get_num_total_reads_for_chromosome(self, chromosome):
        stats = self._get_stats_for_chromosome(chromosome)
        return stats.num_reads

    def get_length_of_chromosome(self, chromosome):
        stats = self._get_stats_for_chromosome(chromosome)
        return stats.chromosome_length


def load_bam_index_stats(bam_file):
    index_stats_text = pysam.idxstats(bam_file.filename)
    lines = index_stats_text.splitlines()
    bam_index_stats = BamIndexStats()

    for line in lines:
        if line.startswith("#"):
            continue

        chromosome, length, num_mapped_reads, num_unmapped_reads = line.split("\t")
        length = int(length)
        num_mapped_reads = int(num_mapped_reads)
        num_unmapped_reads = int(num_unmapped_reads)

        stats_row = BamIndexStatsRow(
            chromosome,
            length,
            num_unmapped_reads + num_mapped_reads,
            num_mapped_reads,
            num_unmapped_reads
        )

        bam_index_stats.add_row(
            stats_row
        )

    return bam_index_stats


def get_num_mapped_reads_covering_chromosome(bam_file, chrom):
    """
    Return the total number of mapped reads covering the specified chromsome
    in the specified BAM file. This is an optimisation, which makes use of the
    index statistics in the BAM index, which record the nunmber of alignments. This is
    an O(1) operation rather than the O(N) operation of looping through all reads in the
    file with read.tid == chromID.
    """
    index_stats = pysam.idxstats(bam_file.filename).splitlines()

    for line in index_stats:
        if not line.startswith("#"):
            the_chrom, length, num_reads, num_unmapped_reads = line.split("\t")
            num_reads = int(num_reads)
            num_unmapped_reads = int(num_unmapped_reads)

            if the_chrom == chrom:
                return num_reads - num_unmapped_reads
    else:
        raise StandardError("Could not find chromosome {} in pysam index stats".format(
            chrom
        ))


def get_num_unmapped_reads_covering_chromosome(bam_file, chrom):
    """
    Return the total number of unmapped reads covering the specified chromsome
    in the specified BAM file. This is an optimisation, which makes use of the
    index statistics in the BAM index, which record the nunmber of alignments. This is
    an O(1) operation rather than the O(N) operation of looping through all reads in the
    file with read.tid == chromID.
    """
    index_stats = pysam.idxstats(bam_file.filename).splitlines()

    for line in index_stats:
        if not line.startswith("#"):
            the_chrom, length, num_reads, num_unmapped_reads = line.split("\t")
            num_reads = int(num_reads)
            num_unmapped_reads = int(num_unmapped_reads)

            if the_chrom == chrom:
                return num_unmapped_reads
    else:
        raise StandardError("Could not find chromosome {} in pysam index stats".format(
            chrom
        ))


def get_total_num_reads_covering_chromosome(bam_file, chrom):
    """
    Return the total number of reads covering the specified chromsome
    in the specified BAM file. This is an optimisation, which makes use of the
    index statistics in the BAM index, which record the nunmber of alignments. This is
    an O(1) operation rather than the O(N) operation of looping through all reads in the
    file with read.tid == chromID.
    """
    index_stats = pysam.idxstats(bam_file.filename).splitlines()

    for line in index_stats:
        if not line.startswith("#"):
            the_chrom, length, num_reads, num_unmapped_reads = line.split("\t")
            num_reads = int(num_reads)
            num_unmapped_reads = int(num_unmapped_reads)

            if the_chrom == chrom:
                return num_reads
    else:
        raise StandardError("Could not find chromosome {} in pysam index stats".format(
            chrom
        ))


def get_valid_chromosome_name(chrom, bam_file):
    """
    Return a chromosome name that matches the naming convention used in the BAM file, i.e.
    add or remove 'chr' as appropriate.
    """
    chr_prefix = bam_file.references[0].startswith('chr')

    if chr_prefix and not chrom.startswith('chr'):
        return 'chr' + chrom

    if not chr_prefix and chrom.startswith('chr'):
        return chrom[3:]

    return chrom

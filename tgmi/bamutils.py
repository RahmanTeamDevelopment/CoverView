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
        return stats.num_unmapped_reads

    def get_num_total_reads_for_chromosome(self, chromosome):
        stats = self._get_stats_for_chromosome(chromosome)
        return stats.num_reads

    def get_length_of_chromosome(self, chromosome):
        stats = self._get_stats_for_chromosome(chromosome)
        return stats.chromosome_length

    def get_total_mapped_reads_in_bam(self):
        return sum(x.num_mapped_reads for x in self.rows)

    def get_total_unmapped_reads_in_bam(self):
        return sum(x.num_unmapped_reads for x in self.rows)

    def get_total_reads_in_bam(self):
        return sum(x.num_reads for x in self.rows)


def extract_bam_index_stats(index_stats_lines):
    bam_index_stats = BamIndexStats()

    for line in index_stats_lines:
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


def load_bam_index_stats_from_file(bam_file):
    index_stats_text = pysam.idxstats(bam_file.filename)
    lines = index_stats_text.splitlines()
    return extract_bam_index_stats(lines)


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

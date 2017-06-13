"""
Assorted utility functions that don't belong anywhere else
"""

from collections import defaultdict
import logging


_logger = logging.getLogger("coverview")


def min_or_nan(data):
    if len(data) == 0:
        return float('NaN')
    else:
        return min(data)


def max_or_nan(data):
    if len(data) == 0:
        return float('NaN')
    else:
        return max(data)


class GenomicInterval(object):
    """
    Represents and interval on the genome. All coordinates are 0-based, and the
    interval is half-open, i.e. the end position is not included in the interval.
    """
    def __init__(self, chromosome, start_pos, end_pos, name=None):
        assert start_pos >= 0
        assert end_pos >= 0

        self.chromosome = chromosome
        self.start_pos = start_pos
        self.end_pos = end_pos

        if name is None:
            self.name = "{}:{}-{}".format(
                chromosome,
                start_pos,
                end_pos
            )
        else:
            self.name = name

    def __eq__(self, other):
        if self.chromosome != other.chromosome:
            return False

        if self.start_pos != other.start_pos:
            return False

        if self.end_pos != other.end_pos:
            return False

        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def size(self):
        return self.end_pos - self.start_pos

    def overlap(self, other):
        if self.chromosome != other.chromosome:
            return GenomicInterval(None, 0, 0)
        else:
            overlap_start = max(self.start_pos, other.start_pos)
            overlap_end = min(self.end_pos, other.end_pos)

            if overlap_end > overlap_start:
                return GenomicInterval(
                    self.chromosome,
                    overlap_start,
                    overlap_end
                )
            else:
                return GenomicInterval(None, 0, 0)


class BedFileParser(object):
    def __init__(self, bed_file):
        self.bed_file = bed_file

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.bed_file.close()

    def next(self):
        line = self.bed_file.next()
        cols = line.strip().split("\t")
        chromosome = cols[0]
        start_pos = int(cols[1])
        end_pos = int(cols[2])
        name = None

        if len(cols) > 3:
            name = cols[4]

        yield GenomicInterval(chromosome, start_pos, end_pos, name)


def get_clusters_of_regions_from_bed_file(bed_file, size_limit=100000):
    """
    Reads a BED file and yields lists of regions that are close
    together.
    """
    all_regions = []

    for line in bed_file:

        line = line.rstrip()

        if len(line) == 0 or line.startswith('#'):
            continue

        row = line.split('\t')

        if len(row) < 4:
            raise StandardError("Error: incorrect BED file!")

        chrom = row[0]
        begin = int(row[1])
        end = int(row[2])
        key = row[3]
        region = "{}:{}-{}".format(chrom, begin, end)
        all_regions.append((chrom, begin, end, region, key))

    if len(all_regions) == 0:
        raise StandardError("No regions in BED file")

    all_regions.sort()
    current_cluster = []

    for chrom, begin, end, region, key in all_regions:
        if len(current_cluster) == 0:
            current_cluster.append((chrom, begin, end, region, key))
        elif current_cluster[-1][0] != chrom:
            yield current_cluster
            current_cluster = [(chrom, begin, end, region, key)]
        elif end - current_cluster[0][1] > size_limit:
            yield current_cluster
            current_cluster = [(chrom, begin, end, region, key)]
        else:
            current_cluster.append((chrom, begin, end, region, key))

    yield current_cluster


def get_names_of_target_regions(bed_file):
    """
    Reads a BED file of target regions and returns a list of target
    names. If there are duplicated names (e.g. multiple exons for the same
    gene with the gene name as the target name) then a number is appended to the
    name e.g. BRCA1_1, BRCA1_2 etc.
    """
    target_names = []
    name_counts = defaultdict(int)

    for line in bed_file:
        line = line.strip()

        if len(line) == 0 or line.startswith("#"):
            continue

        cols = line.split('\t')
        target_name = cols[3]
        name_counts[target_name] += 1

    if len(name_counts) == 0:
        raise StandardError("No regions in BED file")

    for target_name, count in name_counts.iteritems():
        if count == 1:
            target_names.append(target_name)
        else:
            for i in xrange(count):
                target_names.append(
                    "{}_{}".format(target_name, i+1)
                )

    return target_names, len(name_counts)

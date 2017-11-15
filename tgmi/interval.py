"""
Assorted utility functions that don't belong anywhere else
"""


from collections import defaultdict
import functools


@functools.total_ordering
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

    def __lt__(self, other):
        return (self.chromosome, self.start_pos, self.end_pos) < (other.chromosome, other.start_pos, other.end_pos)

    def __hash__(self):
        return hash((
            self.chromosome,
            self.start_pos,
            self.end_pos,
            self.name
        ))

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

    def __repr__(self):
        return "{}\t{}:{}-{}".format(
            self.name,
            self.chromosome,
            self.start_pos,
            self.end_pos
        )

    def __str__(self):
        return self.__repr__()

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


def uniquify_region_names(regions):
    """
    If there are duplicate names in the region list, e.g. 2 regions called BRCA, then
    rename to BRCA_1, BRCA_2 etc. Return a new list of new intervals with unique names.
    """
    regions_by_name = defaultdict(list)
    regions_with_unique_names = []

    for region in sorted(regions):
        regions_by_name[region.name].append(region)

    for name, regions_with_same_name in regions_by_name.iteritems():
        count = len(regions_with_same_name)

        if count == 1:
            old_region = regions_with_same_name[0]

            new_region = GenomicInterval(
                old_region.chromosome,
                old_region.start_pos,
                old_region.end_pos,
                old_region.name
            )

            regions_with_unique_names.append(new_region)
        else:
            for i in xrange(count):
                old_region = regions_with_same_name[i]

                new_region = GenomicInterval(
                    old_region.chromosome,
                    old_region.start_pos,
                    old_region.end_pos,
                    "{}_{}".format(old_region.name, i+1)
                )

                regions_with_unique_names.append(new_region)

    return sorted(regions_with_unique_names)


def cluster_genomic_intervals(intervals, size_limit=100000):
    """
    Reads a BED file and yields lists of regions that are close
    together.
    """
    if len(intervals) == 0:
        raise StandardError("Passed empty list of intervals to cluster function")

    all_intervals = sorted(intervals)
    current_cluster = []

    for interval in all_intervals:
        if len(current_cluster) == 0:
            current_cluster.append(interval)
        elif current_cluster[-1].chromosome != interval.chromosome:
            yield current_cluster
            current_cluster = [interval]
        elif interval.end_pos - current_cluster[0].start_pos > size_limit:
            yield current_cluster
            current_cluster = [interval]
        else:
            current_cluster.append(interval)

    yield current_cluster

"""
Assorted utility functions that don't belong anywhere else
"""

from collections import defaultdict


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


def get_clusters_of_regions_from_bed_file(bed_file_name, size_limit=100000):
    """
    Reads a BED file and yields lists of regions that are close
    together.
    """
    all_regions = []

    with open(bed_file_name, 'r') as bed_file:
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


def get_names_of_target_regions(input_file_name):
    """
    Reads a BED file of target regions and returns a list of target
    names. If there are duplicated names (e.g. multiple exons for the same
    gene with the gene name as the target name) then a number is appended to the
    name e.g. BRCA1_1, BRCA1_2 etc.
    """
    target_names = []
    name_counts = defaultdict(int)

    with open(input_file_name, 'r') as bed_file:
        for line in bed_file:
            line = line.strip()

            if len(line) == 0 or line.startswith("#"):
                continue

            cols = line.split('\t')
            target_name = cols[3]
            name_counts[target_name] += 1

        for target_name, count in name_counts.iteritems():
            if count == 1:
                target_names.append(target_name)
            else:
                for i in xrange(1, count + 1):
                    target_names.append(
                        "{}_{}".format(target_name, i)
                    )

    return target_names, len(name_counts)


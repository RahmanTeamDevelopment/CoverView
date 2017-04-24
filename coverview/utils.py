"""
Assorted utility functions that don't belong anywhere else
"""


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

            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

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


def makeNames(inputf):
    ret = []
    latest = dict()
    for line in open(inputf):
        line = line.rstrip()
        if line == '': continue
        if line.startswith('#'): continue
        cols = line.split('\t')

        if cols[3] in latest.keys():
            latest[cols[3]] += 1
        else:
            latest[cols[3]] = 1

        ret.append(cols[3] + '_' + str(latest[cols[3]]))

    for i in range(len(ret)):
        x = ret[i]
        if x.endswith('_1'):
            if not x[:-2] + '_2' in ret:
                ret[i] = x[:-2]

    return ret, len(latest.keys())


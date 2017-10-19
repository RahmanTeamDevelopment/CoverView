import json

def read_data(prefix):
    """Read all CoverView output data required by the GUI"""

    profiles, region_coords = read_profiles_data(prefix)

    return {
        'regions': read_regions_data(prefix),
        'profiles': profiles,
        'region_coords': region_coords,
        'summary': read_summary_data(prefix)
    }


def read_metadata(prefix):

    with open(prefix + '_meta.json') as json_data:
        ret = json.load(json_data)
        json_data.close()

    ret['config_opts']['count_duplicate_reads'] = 'Included' if ret['config_opts']['count_duplicate_reads'] else 'Excluded'

    return ret


def read_regions_data(prefix):
    """Read _region output file into dict"""

    ret = []
    idx = {}
    columns = []
    for line in open(prefix+'_regions.txt'):
        line = line.strip()
        if line == '':
            continue

        if line.startswith('#'):
            line = line[1:]
            header = line.split()
            for i in range(len(header)):
                idx[header[i]] = i

            columns = ['RC', 'MEDCOV', 'MINCOV', 'MEDQCOV', 'MINQCOV', 'MAXFLMQ', 'MAXFLBQ', 'Pass_or_fail']
            if 'MEDCOV+' in header:
                columns += [
                    'RC+', 'MEDCOV+', 'MINCOV+', 'MEDQCOV+', 'MINQCOV+', 'MAXFLMQ+', 'MAXFLBQ+',
                    'RC-', 'MEDCOV-', 'MINCOV-', 'MEDQCOV-', 'MINQCOV-', 'MAXFLMQ-', 'MAXFLBQ-'
                ]
            continue

        cols = line.split()
        record = {'region': cols[0]}
        for c in columns:
            value = cols[idx[c]]
            if value == '.':
                value = '--'
            elif c.startswith('RC') or c.startswith('MIN'):
                value = int(value)
            elif c != 'Pass_or_fail':
                value = float(value)
            if c.startswith('MED') and value == int(value):
                value = int(value)
            record[c.lower()] = value
        ret.append(record)

    return ret


def read_profiles_data(prefix):
    """Read _profiles output file into dict"""

    ret = {}
    ret_coords = {}
    idx = {}
    columns = []
    region_data = {}
    region = None
    region_chrom = None
    region_start = None
    region_end = None
    for line in open(prefix+'_profiles.txt'):
        line = line.strip()
        if line == '':
            continue

        if line[0] == '#':
            line = line[1:]
            header = line.split()
            for i in range(len(header)):
                idx[header[i]] = i

            columns = ['COV', 'QCOV', 'MEDBQ', 'FLBQ', 'MEDMQ', 'FLMQ']
            if 'COV+' in header:
                columns += [
                    'COV+', 'QCOV+', 'MEDBQ+', 'FLBQ+', 'MEDMQ+', 'FLMQ+',
                    'COV-', 'QCOV-', 'MEDBQ-', 'FLBQ-', 'MEDMQ-', 'FLMQ-'
                ]
            continue

        if line[0] == '[':
            if region is not None:
                ret[region] = region_data
                ret_coords[region] = [region_chrom, int(region_start), int(region_end)+1]
            region = line[1:-1]
            region_data = {}
            region_chrom = region_start = region_end = None
            continue

        cols = line.split()

        if region_start is None:
            region_chrom = cols[0]
            region_start = cols[1]
        region_end = cols[1]

        for c in columns:
            if c not in region_data:
                region_data[c] = []

            value = cols[idx[c]]
            if value == '.':
                value = None
            else:
                value = float(value)

            region_data[c].append(value)

    ret[region] = region_data
    ret_coords[region] = [region_chrom, int(region_start), int(region_end)+1]
    return ret, ret_coords


def read_summary_data(prefix):
    """Read _summary output file into dict"""

    ret = {}
    for line in open(prefix+'_summary.txt'):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        cols = line.split()
        ret[cols[0]] = {'RC': cols[1], 'RCIN': cols[2], 'RCOUT': cols[3]}
    return ret




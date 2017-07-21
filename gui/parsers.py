

def read_data(prefix):
    """Read all CoverView output data required by the GUI"""

    return {
        'regions': read_regions_data(prefix),
        'profiles': read_profiles_data(prefix),
        'summary': read_summary_data(prefix)
    }


def read_regions_data(prefix):
    """Read _region output file into dict"""

    ret = {}
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
                    'MEDCOV+', 'MINCOV+', 'MEDQCOV+', 'MINQCOV+', 'MAXFLMQ+', 'MAXFLBQ+',
                    'MEDCOV-', 'MINCOV-', 'MEDQCOV-', 'MINQCOV-', 'MAXFLMQ-', 'MAXFLBQ-'
                ]
            continue

        cols = line.split()
        record = {}
        for c in columns:
            value = cols[idx[c]]
            if value == '.':
                value = None
            elif value not in ['FAIL', 'PASS']:
                value = float(value)
            record[c] = value
        ret[cols[0]] = record

    return ret


def read_profiles_data(prefix):
    """Read _profiles output file into dict"""

    ret = {}
    idx = {}
    columns = []
    region_data = {}
    region = None
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
            region = line[1:-1]
            region_data = {}
            continue

        cols = line.split()
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
    return ret


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




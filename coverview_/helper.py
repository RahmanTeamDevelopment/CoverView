
def parse_ini_file(fn):
    """Parse INI file and return dictionary"""

    ret = {}
    section = ''
    with open(fn) as infile:
        for line in infile:
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue

            if line[0] == '[' and line[-1] == ']':
                s = line[1:-1]
                if '[' in s or ']' in s or '=' in s:
                    continue
                section = s

            elif line.count('=') == 1:
                key, value = line.split('=')
                key = key.strip()
                value = value.strip()
                if section != '':
                    key = '{}.{}'.format(section, key)
                ret[key.upper()] = value
    return ret


def is_boolean(x):
    return x.upper() in ['TRUE', 'FALSE']


def is_int(x):
    try:
        int(x)
        return True
    except:
        return False


def process_option(_logger, ini_data, key, type, default):
    name = '['+key.lower().replace('.', ']/')
    if key in ini_data:
        if type == 'boolean' and is_boolean(ini_data[key]):
            return ini_data[key].upper() == 'TRUE'
        elif type == 'int' and is_int(ini_data[key]):
            return int(ini_data[key])
        else:
            msg = 'Configuration option \"{}\" has incorrect value ({})'.format(name, ini_data[key])
            _logger.error(msg)
            raise StandardError(msg)
    else:
        return default


def process_pass_option(_logger, ini_data):
    ret = None
    for k, v in ini_data.iteritems():
        if '.' not in k:
            continue
        [section, flag] = k.split('.')
        if section != 'PASS':
            continue
        if flag not in ['MINCOV_MIN', 'MINQCOV_MIN', 'MAXFLBQ_MAX', 'MAXFLMQ_MAX', 'MEDCOV_MIN', 'MEDQCOV_MIN']:
            continue
        if ret is None:
            ret = {}
        try:
            ret[flag] = int(v)
        except:
            msg = 'Configuration option \"[pass]/{}\" has incorrect value ({})'.format(flag, ini_data['PASS.'+flag])
            _logger.error(msg)
            raise StandardError(msg)
    return ret


def read_config_file(fn, _logger):
    ret = {}

    ini_data = parse_ini_file(fn)

    ret['count_duplicate_reads'] = process_option(_logger, ini_data, 'READS.DUPLICATES', 'boolean', True)
    ret['direction'] = process_option(_logger, ini_data, 'READS.DIRECTION', 'boolean', False)
    ret['only_flagged_profiles'] = process_option(_logger, ini_data, 'OUTPUTS.ONLY_FLAGGED_PROFILES', 'boolean', False)
    ret['low_bq'] = process_option(_logger, ini_data, 'QUALITY.LOW_BQ', 'int', 10)
    ret['low_mq'] = process_option(_logger, ini_data, 'QUALITY.LOW_MQ', 'int', 20)

    ret['outputs'] = {}
    ret['outputs']['regions'] = process_option(_logger, ini_data, 'OUTPUTS.REGIONS_FILE', 'boolean', True)
    ret['outputs']['profiles'] = process_option(_logger, ini_data, 'OUTPUTS.PROFILES_FILE', 'boolean', True)

    ret['transcript'] = {}
    ret['transcript']['regions'] = process_option(_logger, ini_data, 'TRANSCRIPT.REGIONS_FILE', 'boolean', True)
    ret['transcript']['profiles'] = process_option(_logger, ini_data, 'TRANSCRIPT.PROFILES_FILE', 'boolean', True)

    ret['pass'] = process_pass_option(_logger, ini_data)

    return ret
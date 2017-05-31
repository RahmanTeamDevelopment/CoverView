import csv


def float_or_nan(x):
    if x == ".":
        return float('NaN')
    else:
        return float(x)


_value_type_map = {
    "#Region": str,
    "Chromosome": lambda x: str(x).replace("chr", "").replace("CHR", ""),
    "Start_position": int,
    "End_position": int,
    "Pass_or_fail": str,
    "RC": int,
    "MEDCOV": float_or_nan,
    "MINCOV": float_or_nan,
    "MEDQCOV": float_or_nan,
    "MINQCOV": float_or_nan,
    "MAXFLMQ": float_or_nan,
    "MAXFLBQ": float_or_nan
}


def load_coveriew_regions_output(file_name):

    data = {}

    with open(file_name, 'r') as regions_file:
        regions_file_reader = csv.DictReader(regions_file, delimiter='\t')

        for record in regions_file_reader:
            region_name = record['#Region']
            data[region_name] = {}

            for key, converter_function in _value_type_map.iteritems():

                if key in record:
                    typed_value = converter_function(record[key])
                    data[region_name][key] = typed_value

    return data

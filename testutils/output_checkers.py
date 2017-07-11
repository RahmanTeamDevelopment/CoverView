import csv


def float_or_nan(x):
    if x == ".":
        return x
    else:
        return float(x)


_regions_value_type_map = {
    "#Region": str,
    "Chromosome": lambda x: str(x).replace("chr", "").replace("CHR", ""),
    "Start_position": int,
    "End_position": int,
    "Start_transcript": str,
    "End_transcript": str,
    "Pass_or_fail": str,
    "RC": int,
    "MEDCOV": float_or_nan,
    "MINCOV": float_or_nan,
    "MEDQCOV": float_or_nan,
    "MINQCOV": float_or_nan,
    "MAXFLMQ": float_or_nan,
    "MAXFLBQ": float_or_nan
}


_summary_value_type_map = {
    "#CHROM": lambda x: x,
    "RC": int,
    "RCIN": lambda x: x if x == "-" else int(x),
    "RCOUT": lambda x: x if x == "-" else int(x)
}


_profile_value_type_map = {
    "#Chromosome": lambda x: x,
    "Position": int,
    "COV": int,
    "QCOV": int,
    "MEDBQ": float_or_nan,
    "FLBQ": float_or_nan,
    "MEDMQ": float_or_nan,
    "FLMQ": float_or_nan
}


def load_coverview_profile_output(file_name):

    data = {}
    header_names = []
    current_region_name = None

    with open(file_name, 'r') as profile_file:
        for line in profile_file:
            if line.startswith("#"):
                header_names = line.strip().split("\t")
            elif line.startswith("["):
                region_name = line.strip().strip("[").strip("]")
                data[region_name] = {}
                current_region_name = region_name
            else:

                if line.strip() == "":
                    continue

                cols = line.strip().split("\t")
                chrom_pos = ":".join([cols[0], cols[1]])
                data[current_region_name][chrom_pos] = {}

                for col_index, col_data in enumerate(cols):
                    typed_data = _profile_value_type_map[header_names[col_index]](col_data)
                    data[current_region_name][chrom_pos][header_names[col_index]] = typed_data
    return data


def load_coverview_summary_output(file_name):
    with open(file_name, 'r') as summary_file:
        summary_file_reader = csv.DictReader(summary_file, delimiter='\t')
        data = {}

        for record in summary_file_reader:
            chrom = record['#CHROM']
            data[chrom] = {}

            for key, value in record.items():
                if key == "#CHROM":
                    continue
                else:
                    typed_value = _summary_value_type_map[key](value)
                    data[chrom][key] = typed_value

        return data


def load_coveriew_regions_output(file_name):
    with open(file_name, 'r') as regions_file:
        data = {}
        regions_file_reader = csv.DictReader(regions_file, delimiter='\t')

        for record in regions_file_reader:
            region_name = record['#Region']
            data[region_name] = {}

            for key, converter_function in _regions_value_type_map.items():

                if key in record:
                    typed_value = converter_function(record[key])
                    data[region_name][key] = typed_value

        return data

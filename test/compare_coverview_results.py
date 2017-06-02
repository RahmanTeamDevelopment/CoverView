#!env/bin/python

"""
Compares the results of the CoverView regions output with a test file. Every value in the
2 files will be compared, and a digest of the differences printed to the command line
"""

import argparse
import csv
import logging


_logger = logging.getLogger("coverview_test")

_value_type_map = {
    "#Region": str,
    "Chromosome": lambda x: str(x).replace("chr", "").replace("CHR", ""),
    "Start_position": int,
    "End_position": int,
    "Pass_or_fail": str,
    "RC": int,
    "MEDCOV": float,
    "MINCOV": float,
    "MEDQCOV": float,
    "MINQCOV": float,
    "MAXFLMQ": float,
    "MAXFLBQ": float
}


def configure_logging():
    """
    Currently just logging to the terminal stderr stream, but this could easily be extended
    to produce a log file or e.g. email alerts.
    """
    logger = logging.getLogger("coverview_test")

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(pathname)s - Line %(lineno)s - %(message)s")

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)

    logger.addHandler(stream_handler)
    logger.setLevel(logging.DEBUG)


def parse_command_line_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--new_file",
        dest='new_file',
        action='store',
        help="File containing new 'regions' output, to be checked against the old output",
        required=True
    )

    parser.add_argument(
        "--old_file",
        dest='old_file',
        action='store',
        help="File containing old 'regions' output. This will be used as the benchmark for testing",
        required=True
    )

    options = parser.parse_args()
    return options


def compare_regions_file_contents(old_file, new_file):
    old_data_by_region = {}

    for old_record in old_file:
        old_data_by_region[old_record['#Region']] = old_record

    for new_record in new_file:
        region_name = new_record['#Region']

        if region_name not in old_data_by_region:
            _logger.error("Region {} does not exist in the test dataset but is in the new dataset".format(
                region_name
            ))
            continue

        old_record = old_data_by_region[region_name]

        for key, old_value in old_record.iteritems():
            if key not in new_record:
                _logger.error("Data value {} is missing from record {}".format(
                    key,
                    region_name
                ))
            else:
                new_value = new_record[key]

                typed_old_value = _value_type_map[key](old_value)
                typed_new_value = _value_type_map[key](new_value)

                if typed_old_value != typed_new_value:

                    if _value_type_map[key] != str:

                        if (abs(typed_old_value - typed_new_value) / typed_old_value) > 0.05:
                            _logger.error("New value for data point {}:{} ({}) differs by more than 1% from old value ({})".format(
                                region_name,
                                key,
                                new_value,
                                old_value
                            ))
                    else:
                        _logger.error(
                            "New value for data point {}:{} ({}) != old value ({})".format(
                                region_name,
                                key,
                                new_value,
                                old_value
                            ))


if __name__ == "__main__":
    configure_logging()
    _logger.info("Parsing command-line arguments")
    options = parse_command_line_arguments()
    _logger.info(options)

    with open(options.old_file, 'r') as old_file:
        with open(options.new_file, 'r') as new_file:
            old_file_reader = csv.DictReader(old_file, delimiter='\t')
            new_file_reader = csv.DictReader(new_file, delimiter='\t')

            compare_regions_file_contents(old_file_reader, new_file_reader)
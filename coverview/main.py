from __future__ import division

import argparse
import collections
import json
import logging
import pysam
import tgmi.bed
import tgmi.interval
import tgmi.math

from . import output
from .calculators import calculate_chromosome_coverage_metrics, get_region_coverage_summary
from .calculators import calculate_minimal_chromosome_coverage_metrics
from .statistics import median


_version = 'v1.3.0'
_logger = logging.getLogger("coverview")


class CoverageCalculator(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.bam_file = pysam.Samfile(options.input, "rb")
        self.transcript_database = None
        self.out_poor = None
        self.num_reads_on_target = collections.defaultdict(int)
        self.ids_of_failed_targets = set()
        self.regions_output = None
        self.per_base_output = None

        if options.transcript_db is not None:
            self.transcript_database = pysam.Tabixfile(
                options.transcript_db
            )

        self.first = True

        if config['outputs']['regions']:
            self.regions_output = output.RegionsOutput(options, config)

        if config['outputs']['profiles']:
            self.per_base_output = output.PerBaseCoverageOutput(options, config)

    def does_region_pass_coverage_thresholds(self, target):
        """
        Returns true if the specified region passes all criteria. The thresholds
        are specified in a configuration file as a dictionary where the keys are
        a string of metric name and either MIN or MAX separated by an underscore, e.g.
        MEDQCOV_MIN and the values are the thresholds.
        """
        if self.config['pass'] is None:
            return True
        else:
            region_coverage_summary = target.summary

            for key, value in self.config['pass'].items():
                metric_name, min_or_max = key.split('_')
                metric_value = region_coverage_summary[metric_name]
                threshold = float(value)

                if min_or_max == 'MIN' and metric_value < threshold:
                    return False

                if min_or_max == 'MAX' and metric_value > threshold:
                    return False

            return True

    def compute_summaries_of_region_coverage(self, profile):
        """
        Caclculates medians, minimums and maximums of various coverage metrics
        for a single region.
        """
        summary = {}

        summary['MEDCOV'] = median(profile.coverage_at_each_base)
        summary['MEDQCOV'] = median(profile.high_quality_coverage_at_each_base)
        summary['MINCOV'] = int(tgmi.math.min_or_nan(profile.coverage_at_each_base))
        summary['MINQCOV'] = int(tgmi.math.min_or_nan(profile.high_quality_coverage_at_each_base))
        summary['MAXFLBQ'] = round(tgmi.math.max_or_nan(profile.fraction_of_low_base_qualities_at_each_base), 3)
        summary['MAXFLMQ'] = round(tgmi.math.max_or_nan(profile.fraction_of_low_mapping_qualities_at_each_base), 3)

        if self.config['direction']:
            summary['MEDCOV_f'] = median(profile.forward_coverage_at_each_base)
            summary['MEDQCOV_r'] = median(profile.reverse_high_quality_coverage_at_each_base)
            summary['MEDQCOV_f'] = median(profile.forward_high_quality_coverage_at_each_base)
            summary['MEDCOV_r'] = median(profile.reverse_coverage_at_each_base)
            summary['MINCOV_f'] = int(tgmi.math.min_or_nan(profile.forward_coverage_at_each_base))
            summary['MINQCOV_f'] = int(tgmi.math.min_or_nan(profile.forward_high_quality_coverage_at_each_base))

            summary['MAXFLBQ_f'] = round(
                tgmi.math.max_or_nan(profile.forward_fraction_of_low_base_qualities_at_each_base),
                3
            )
            summary['MAXFLMQ_f'] = round(
                tgmi.math.max_or_nan(profile.forward_fraction_of_low_mapping_qualities_at_each_base),
                3
            )

            summary['MINCOV_r'] = int(tgmi.math.min_or_nan(profile.reverse_coverage_at_each_base))
            summary['MINQCOV_r'] = int(tgmi.math.min_or_nan(profile.reverse_high_quality_coverage_at_each_base))
            summary['MAXFLBQ_r'] = round(
                tgmi.math.max_or_nan(profile.reverse_fraction_of_low_base_qualities_at_each_base),
                3
            )

            summary['MAXFLMQ_r'] = round(
                tgmi.math.max_or_nan(profile.reverse_fraction_of_low_mapping_qualities_at_each_base),
                3
            )

        return summary

    def close_output_files(self):
        """
        Make sure that all output files are closed. The way CoverView is configured means that there
        are several optional output files, which may or may not exist depending on flags in the config
        file.
        """
        if self.config['outputs']['profiles']:
            self.out_profiles.close()
            if self.transcript_database is not None:
                self.out_poor.close()

    def write_output_file_headers(self):
        """
        Write the header rows for each output file, if required.
        """
        _logger.info("Writing output headers")

        if self.regions_output is not None:
            self.regions_output.write_header()

        if self.per_base_output is not None:
            self.per_base_output.write_header()

    def write_outputs_for_region(self, region_coverage_data):
        """
        Write various summaries for each targeted region, if required.
        """
        if self.regions_output:
            self.regions_output.write_output(
                region_coverage_data,
                self.transcript_database
            )

        if self.per_base_output:
            self.per_base_output.write_output(
                region_coverage_data,
                self.transcript_database
            )

    def calculate_coverage_summaries(self, intervals):
        _logger.info("Coverage metrics will be generated in a single process")
        self.write_output_file_headers()
        num_clusters = 0

        for cluster in tgmi.interval.cluster_genomic_intervals(intervals):
            num_clusters += 1
            for target in get_region_coverage_summary(self.bam_file, cluster, self.config):

                if target is None:
                    continue

                per_base_summary = target.per_base_coverage_profile
                self.num_reads_on_target[target.chromosome] += per_base_summary.num_reads_in_region

                target.summary = self.compute_summaries_of_region_coverage(
                    target.per_base_coverage_profile
                )

                target.passes_thresholds = self.does_region_pass_coverage_thresholds(
                    target
                )

                if not target.passes_thresholds:
                    if '_' in target.region_name:
                        ids = target.region_name[:target.region_name.find('_')]
                    else:
                        ids = target.region_name

                    self.ids_of_failed_targets.add(ids)
                self.write_outputs_for_region(target)

        _logger.info("Finished computing coverage metrics in all regions")
        _logger.debug("Data was processed in {} clusters".format(num_clusters))


def get_default_config():
    return {
        "count_duplicate_reads": True,
        "outputs": {
            "regions": True,
            "profiles": True,
        },
        "transcript": {
            "regions": True,
            "profiles": True
        },
        "low_bq": 10,
        "low_mq": 20,
        "only_fail_profiles": False,
        "pass": None,
        "direction": False,
    }


def load_and_validate_config(config_file_name):
    """
    Part of the command line input is a configuration file in JSON format. Here
    we validate the contents of the JSON file.
    """
    config = get_default_config()
    input_config = None

    allowed_config_parameters = {
        "count_duplicate_reads",
        "direction",
        "low_bq",
        "low_mq",
        "only_fail_profiles",
        "outputs",
        "pass",
        "transcript",
    }

    allowed_config_paramters_outputs = {
        "profiles",
        "regions",
        "summary"
    }

    if config_file_name is not None:
        with open(config_file_name) as config_file:

            try:
                input_config = json.load(config_file)
                _logger.debug(input_config)
            except:
                _logger.error("Invalid JSON config file")
                _logger.error("File {} cannot be loaded with the Pyton JSON parser".format(config_file_name))
                _logger.error("Check the file for JSON format errors")

        for key, value in input_config.items():

            if key not in allowed_config_parameters:
                _logger.error("Invalid parameter '{}' found in config JSON file".format(key))
                raise StandardError("Invalid configuration file")

            if isinstance(value, collections.Mapping):
                if key not in config or config[key] is None:
                    config[key] = {}

                for key_2, value_2 in value.items():

                    if key == "outputs" and key_2 not in allowed_config_paramters_outputs:
                        _logger.error("Invalid outputs parameter '{}' found in config JSON file".format(key_2))
                        raise StandardError("Invalid outputs section in configuration file")
                    else:
                        config[key][key_2] = value_2
            else:
                config[key] = value

    return config


def get_input_options(command_line_args):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        dest='input',
        action='store',
        help="Input (BAM) filename",
        required=True
    )

    parser.add_argument(
        "-o",
        "--output",
        default='output',
        dest='output',
        action='store',
        help="Output filename"
    )

    parser.add_argument(
        "-b",
        "--bed",
        default=None,
        dest='bedfile',
        action='store',
        help="Input BED filename"
    )

    parser.add_argument(
        "-c",
        "--config",
        default=None,
        dest='config',
        action='store',
        help="Configuration file"
    )

    parser.add_argument(
        "-t",
        "--transcript_db",
        default=None,
        dest='transcript_db',
        action='store',
        help="Tabix-indexed database of transcripts"
    )

    options = parser.parse_args(command_line_args)
    config = load_and_validate_config(options.config)

    return options, config


def configure_logging():
    """
    Currently just logging to the terminal stderr stream, but this could easily be extended
    to produce a log file or e.g. email alerts.
    """
    logger = logging.getLogger("coverview")

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(filename)s - Line %(lineno)s - %(message)s")

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)

    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO)

    logger.info('CoverView {} started running'.format(_version))


def main(command_line_args):
    configure_logging()
    options, config = get_input_options(command_line_args)

    _logger.debug("Running CoverView {} with options".format(_version))
    _logger.debug(options)
    _logger.debug(config)

    bam_file = pysam.Samfile(options.input, "rb")

    if options.bedfile is None:
        _logger.info("No input BED file specified. Computing minimal coverage information")

        chromosome_coverage_metrics = calculate_minimal_chromosome_coverage_metrics(
            bam_file,
            options
        )

        output.output_minimal_chromosome_coverage_metrics(
            options,
            chromosome_coverage_metrics
        )

        _logger.info('CoverView {} succesfully finished'.format(_version))
    else:
        with open(options.bedfile) as bed_file:
            bed_parser = tgmi.bed.BedFileParser(bed_file)
            all_regions = []

            for region in bed_parser:
                all_regions.append(region)

            regions_with_unique_names = tgmi.interval.uniquify_region_names(all_regions)

        number_of_targets = len(regions_with_unique_names)

        _logger.info("There are {} target regions".format(number_of_targets))

        coverage_calculator = CoverageCalculator(options, config)
        coverage_calculator.calculate_coverage_summaries(
            regions_with_unique_names
        )

        ids_of_failed_targets = coverage_calculator.ids_of_failed_targets
        num_failed_targets = len(ids_of_failed_targets)

        _logger.info("{} regions failed the coverage thresholds".format(num_failed_targets))

        chromosome_coverage_metrics = calculate_chromosome_coverage_metrics(
            bam_file,
            coverage_calculator.num_reads_on_target
        )

        output.output_chromosome_coverage_metrics(
            options,
            chromosome_coverage_metrics
        )

        _logger.info("CoverView {} succesfully finished".format(_version))

    return 0  # Standard success code

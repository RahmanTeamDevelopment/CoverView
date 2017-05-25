from __future__ import division

import argparse
import bamgen
import json
import logging
import os
import pysam
import shutil

import warnings

from . import output
from .calculators import calculate_chromosome_coverage_metrics, get_region_coverage_summary
from .calculators import calculate_minimal_chromosome_coverage_metrics
from .statistics import median
from .utils import *


_version = 'v1.2.0'
_logger = logging.getLogger("coverview")


class CoverageCalculator(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.bam_file = pysam.Samfile(options.input, "rb")
        self.transcript_database = None
        self.out_poor = None
        self.num_reads_on_target = {}
        self.ids_of_failed_targets = set()
        self.regions_output = None
        self.per_base_output = None
        self.gui_output = None

        if config['reference_file'] is not None:
            ref_file_name = config['reference_file']

            if ref_file_name == "__MOCK__":
                _logger.info("Using Mock reference file. This should only be used for testing")
                self.ref_file = bamgen.MockReferenceFile()
            else:
                _logger.info("Using file {} to retrieve reference sequence".format(
                    ref_file_name
                ))
                self.ref_file = pysam.Fastafile(ref_file_name)

        if config['transcript_db'] is not None:
            self.transcript_database = pysam.Tabixfile(config['transcript_db'])

        self.first = True

        if config['outputs']['regions']:
            self.regions_output = output.RegionsOutput(options, config)

        if config['outputs']['profiles']:
            self.per_base_output = output.PerBaseCoverageOutput(options, config)

        if config['outputs']['gui']:
            self.gui_output = output.GuiOutput(options, config)

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

            for key, value in self.config['pass'].iteritems():
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
        summary['MINCOV'] = int(min_or_nan(profile.coverage_at_each_base))
        summary['MINQCOV'] = int(min_or_nan(profile.high_quality_coverage_at_each_base))
        summary['MAXFLBQ'] = round(max_or_nan(profile.fraction_of_low_base_qualities_at_each_base), 3)
        summary['MAXFLMQ'] = round(max_or_nan(profile.fraction_of_low_mapping_qualities_at_each_base), 3)

        if self.config['direction']:
            summary['MEDCOV_f'] = median(profile.forward_coverage_at_each_base)
            summary['MEDQCOV_r'] = median(profile.reverse_high_quality_coverage_at_each_base)
            summary['MEDQCOV_f'] = median(profile.forward_high_quality_coverage_at_each_base)
            summary['MEDCOV_r'] = median(profile.reverse_coverage_at_each_base)
            summary['MINCOV_f'] = int(min_or_nan(profile.forward_coverage_at_each_base))
            summary['MINQCOV_f'] = int(min_or_nan(profile.forward_high_quality_coverage_at_each_base))
            summary['MAXFLBQ_f'] = round(max_or_nan(profile.forward_fraction_of_low_base_qualities_at_each_base), 3)
            summary['MAXFLMQ_f'] = round(max_or_nan(profile.forward_fraction_of_low_mapping_qualities_at_each_base), 3)
            summary['MINCOV_r'] = int(min_or_nan(profile.reverse_coverage_at_each_base))
            summary['MINQCOV_r'] = int(min_or_nan(profile.reverse_high_quality_coverage_at_each_base))
            summary['MAXFLBQ_r'] = round(max_or_nan(profile.reverse_fraction_of_low_base_qualities_at_each_base), 3)
            summary['MAXFLMQ_r'] = round(max_or_nan(profile.reverse_fraction_of_low_mapping_qualities_at_each_base), 3)

        return summary

    def close_output_files(self):
        """
        Make sure that all output files are closed. The way CoverView is configured means that there
        are several optional output files, which may or may not exist depending on flags in the config
        file.
        """
        if self.config['outputs']['profiles']:
            self.out_profiles.close()
            if not config['transcript_db'] is None:
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

        if self.gui_output is not None:
            self.gui_output.write_header()

    def write_outputs_for_region(self, region_coverage_data):
        """
        Write various summaries for each targeted region, if required.
        """
        if self.regions_output:
            self.regions_output.write_output(region_coverage_data)

        if self.per_base_output:
            self.per_base_output.write_output(region_coverage_data)

        if self.gui_output:
            self.gui_output.write_output(region_coverage_data)

    def calculate_coverage_summaries(self):
        _logger.info("Coverage metrics will be generated in a single process")
        self.write_output_file_headers()
        num_clusters = 0

        for cluster in get_clusters_of_regions_from_bed_file(self.options.bedfile):
            num_clusters += 1
            for target in get_region_coverage_summary(self.bam_file, cluster, self.config):

                if target is None:
                    continue

                region_name = target.region_name
                per_base_summary = target.per_base_coverage_profile
                self.num_reads_on_target[region_name] = per_base_summary.num_reads_in_region

                target.summary = self.compute_summaries_of_region_coverage(
                    target.per_base_coverage_profile
                )

                target.passes_thresholds = self.does_region_pass_coverage_thresholds(
                    target
                )

                if self.gui_output is not None:
                    target.Ref = self.getReferenceSequence(
                        target.chromosome,
                        target.start_position,
                        target.end_position
                    )
                else:
                    target.Ref = None

                if not target.passes_thresholds:
                    if '_' in target.region_name:
                        ids = target.region_name[:target.region_name.find('_')]
                    else:
                        ids = target.region_name

                    self.ids_of_failed_targets.add(ids)
                self.write_outputs_for_region(target)

        _logger.info("Finished computing coverage metrics in all regions")
        _logger.info("Data was processed in {} clusters".format(num_clusters))

    def getReferenceSequence(self, chrom, start, end):
        start = max(1, start)
        end = min(end, self.ref_file.getReferenceLength(chrom))
        seq = self.ref_file.fetch(chrom, start - 1, end)

        return seq.upper()


def get_default_config():
    return {
        "duplicates": True,
        "outputs": {
            "regions": True,
            "profiles": True,
            "gui": False,
        },
        "transcript": {
            "regions": True,
            "profiles": True
        },
        "low_bq": 10,
        "low_mq": 20,
        "only_fail_profiles": False,
        "transcript_db": None,
        "pass": None,
        "direction": False
    }


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

    options = parser.parse_args(command_line_args)
    config = get_default_config()

    if options.config is not None:
        with open(options.config) as config_file:
            input_config = json.load(config_file)
            for key, value in input_config.iteritems():
                if key in ('outputs', 'transcripts'):
                    for key_2, value_2 in value.iteritems():
                        config[key][key_2] = value_2
                else:
                    config[key] = value

    return options, config


def configure_logging():
    """
    Currently just logging to the terminal stderr stream, but this could easily be extended
    to produce a log file or e.g. email alerts.
    """
    logger = logging.getLogger("coverview")

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(pathname)s - Line %(lineno)s - %(message)s")

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)

    logger.addHandler(stream_handler)
    logger.setLevel(logging.DEBUG)

    warnings.simplefilter("ignore", RuntimeWarning)

    logger.info('CoverView v1.2.0 started running')


def clean_up_old_gui_output(gui_output_html_file, gui_data_directory):
    if os.path.exists(gui_output_html_file):
        os.remove(gui_output_html_file)

    if os.path.exists(gui_data_directory):
        os.removedirs(gui_data_directory)


def create_gui_output_directory(options, config):

    _logger.info("Creating output directory structure for GUI in directory".format(
        config['outputs']['gui_output_directory']
    ))

    template_gui_html_file = config['gui']['template_gui_html_file']
    javascript_directory = config['gui']['javascript_directory']
    gui_output_direcory = config['outputs']['gui_output_directory']
    gui_data_directory = os.path.join(gui_output_direcory, "data")
    gui_output_html_file = os.path.join(gui_output_direcory, "_coverview.html")
    gui_output_javascript_directory = os.path.join(gui_output_direcory, "lib")

    clean_up_old_gui_output(gui_output_html_file, gui_data_directory)
    os.makedirs(os.path.join(gui_output_direcory, 'data'))
    shutil.copy(template_gui_html_file, gui_output_html_file)
    shutil.copytree(javascript_directory, gui_output_javascript_directory)


def main(command_line_args):
    configure_logging()
    options, config = get_input_options(command_line_args)

    _logger.info("Running CoverView {} with options".format(_version))
    _logger.info(command_line_args)
    _logger.info(options)
    _logger.info(config)

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
        if config['outputs']['gui']:
            create_gui_output_directory(options, config)

        target_names, unique_target_ids = get_names_of_target_regions(options.bedfile)
        number_of_targets = len(target_names)

        _logger.info("There are {} target regions".format(number_of_targets))

        coverage_calculator = CoverageCalculator(options, config)
        coverage_calculator.calculate_coverage_summaries()

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

        if config['outputs']['gui']:
            coverage_calculator.gui_output.finalize_output(
                options,
                chromosome_coverage_metrics,
                config,
                number_of_targets,
                num_failed_targets,
                unique_target_ids,
                ids_of_failed_targets
            )

        _logger.info("CoverView {} succesfully finished".format(_version))
        return 0  # Standard success code

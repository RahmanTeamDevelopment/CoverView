#!env/bin/python

from __future__ import division

from optparse import OptionParser

import json
import logging
import numpy
import os
import pysam
import shutil
import warnings

from coverview import output
from coverview.calculators import calculate_chromosome_coverage_metrics, get_region_coverage_summary
from coverview.calculators import calculate_minimal_chromosome_coverage_metrics
from coverview.utils import *

__version = 'v1.2.0'
_logger = logging.getLogger("coverview")


class CoverageCalculator(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.bam_file = pysam.Samfile(options.input, "rb")
        self.transcript_database = None
        self.out_targets = None
        self.out_poor = None
        self.out_json = None
        self.out_profiles = None
        self.num_reads_on_target = {}
        self.ids_of_failed_targets = set()

        if not config['transcript_db'] is None:
            self.transcript_database = pysam.Tabixfile(config['transcript_db'])

        self.first = True

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
            region_coverage_summary = target['Summary']

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

        summary['MEDCOV'] = numpy.median(profile['COV'])
        summary['MEDQCOV'] = numpy.median(profile['QCOV'])
        summary['MINCOV'] = min_or_nan(profile['COV'])
        summary['MINQCOV'] = min_or_nan(profile['QCOV'])
        summary['MAXFLBQ'] = max_or_nan(profile['FLBQ'])
        summary['MAXFLMQ'] = max_or_nan(profile['FLMQ'])

        if self.config['direction']:
            summary['MEDCOV_f'] = numpy.median(profile['COV_f'])
            summary['MEDQCOV_r'] = numpy.median(profile['QCOV_r'])
            summary['MEDQCOV_f'] = numpy.median(profile['QCOV_f'])
            summary['MEDCOV_r'] = numpy.median(profile['COV_r'])
            summary['MINCOV_f'] = min_or_nan(profile['COV_f'])
            summary['MINQCOV_f'] = min_or_nan(profile['QCOV_f'])
            summary['MAXFLBQ_f'] = max_or_nan(profile['FLBQ_f'])
            summary['MAXFLMQ_f'] = max_or_nan(profile['FLMQ_f'])
            summary['MINCOV_r'] = min_or_nan(profile['COV_r'])
            summary['MINQCOV_r'] = min_or_nan(profile['QCOV_r'])
            summary['MAXFLBQ_r'] = max_or_nan(profile['FLBQ_r'])
            summary['MAXFLMQ_r'] = max_or_nan(profile['FLMQ_r'])

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

    def calculate_coverage_summaries(self):
        _logger.info("Coverage metrics will be generated in a single process")
        _logger.info("Writing output headers")

        output.output_target_file_header(
            self.config,
            self.out_poor,
            self.out_json,
            self.out_targets,
            self.out_profiles
        )

        num_clusters = 0

        for cluster in get_clusters_of_regions_from_bed_file(options.bedfile):
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
                    per_base_summary
                )

                if self.config['outputs']['gui']:
                    target.Ref = self.getReferenceSequence(
                        target.chromosome,
                        target.start_position,
                        target.end_position
                    )
                else:
                    coverage_data['Ref'] = None

                for o in self.outputs:
                    o.compute_coverage_metric(target)

                if not target['PASS']:
                    if '_' in target['Name']:
                        ids = target['Name'][:target['Name'].find('_')]
                    else:
                        ids = target['Name']

                    self.ids_of_failed_targets.append(ids)
                self.output(target)
        self.close_output_files()

        _logger.info("Finished computing coverage metrics in all regions")
        _logger.info("Data was processed in {} clusters".format(num_clusters))

    def getReferenceSequence(self, chrom, start, end):
        start = max(1, start)
        end = min(end, self.reffile.getReferenceLength(chrom))
        seq = self.reffile.fetch(chrom, start - 1, end)

        return seq.upper()

    def output(self, target):
        if self.config['outputs']['regions']:
            self.output_target(target)

        if self.config['outputs']['profiles']:
            if not self.config['only_fail_profiles']:
                output.output_profiles(
                    target,
                    self.out_profiles
                )
            else:
                if not self.config['pass'] is None:
                    if not target['PASS']:
                        coverage.output.output_profiles(
                            target,
                            self.out_profiles
                        )
                else:
                    coverage.output.output_profiles(
                        target,
                        self.out_profiles
                    )

        if self.config['outputs']['gui']:
            self.output_json(target)

        if self.first:
            self.first = False

    def output_json(self, target):
        if not self.first:
            self.out_json.write(',')

        self.out_json.write(json.dumps(target, separators=(',', ':')))

    def output_target(self, target):

        summary = target['Summary']

        record = [
            target['Name'],
            target['Chrom'],
            str(target['Start']),
            str(target['End'])
        ]

        if self.config['transcript']['regions'] and not self.config['transcript_db'] is None:

            chrom = target['Chrom']
            region_start = target['Start']
            region_end = target['End']

            transcripts_overlapping_start = get_transcripts_overlapping_position(
                self.transcript_database, chrom, region_start
            )

            transcripts_overlapping_end = get_transcripts_overlapping_position(
                self.transcript_database, chrom, region_end
            )

            record.extend([
                transcripts_overlapping_start,
                transcripts_overlapping_end
            ])

        if not self.config['pass'] is None:
            if target['PASS']:
                record.append('PASS')
            else:
                record.append('FAIL')

        if str(summary['MAXFLMQ']) == 'nan':
            summary['MAXFLMQ'] = '.'

        if str(summary['MAXFLBQ']) == 'nan':
            summary['MAXFLBQ'] = '.'

        record.extend([
            str(target['Profiles']['RC']),
            str(summary['MEDCOV']),
            str(summary['MINCOV']),
            str(summary['MEDQCOV']),
            str(summary['MINQCOV']),
            str(summary['MAXFLMQ']),
            str(summary['MAXFLBQ'])
        ])

        if self.config['direction']:
            if str(summary['MAXFLMQ_f']) == 'nan':
                summary['MAXFLMQ_f'] = '.'

            if str(summary['MAXFLBQ_f']) == 'nan':
                summary['MAXFLBQ_f'] = '.'

            if str(summary['MAXFLMQ_r']) == 'nan':
                summary['MAXFLMQ_r'] = '.'

            if str(summary['MAXFLBQ_r']) == 'nan':
                summary['MAXFLBQ_r'] = '.'

            record.extend([
                str(summary['MEDCOV_f']),
                str(summary['MINCOV_f']),
                str(summary['MEDQCOV_f']),
                str(summary['MINQCOV_f']),
                str(summary['MAXFLMQ_f']),
                str(summary['MAXFLBQ_f'])
            ])

            record.extend([
                str(summary['MEDCOV_r']),
                str(summary['MINCOV_r']),
                str(summary['MEDQCOV_r']),
                str(summary['MINQCOV_r']),
                str(summary['MAXFLMQ_r']),
                str(summary['MAXFLBQ_r'])
            ])

        self.out_targets.write('\t'.join(record) + '\n')


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


def get_input_options():
    parser = OptionParser(usage='usage: python %prog [options]')

    parser.add_option(
        "-i",
        "--input",
        default='input.bam',
        dest='input',
        action='store',
        help="Input (BAM) filename [default value: %default]"
    )

    parser.add_option(
        "-o",
        "--output",
        default='output',
        dest='output',
        action='store',
        help="Output filename [default value: %default]"
    )

    parser.add_option(
        "-b",
        "--bed",
        default=None,
        dest='bedfile',
        action='store',
        help="Input BED filename [default value: %default]"
    )

    parser.add_option(
        "-c",
        "--config",
        default=None,
        dest='config',
        action='store',
        help="Configuration file [default value: %default]"
    )

    options, args = parser.parse_args()
    config = get_default_config()

    if options.config is not None:
        with open(options.config) as config_file:
            input_config = json.load(config_file)
            for key,value in input_config.iteritems():
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

    # TODO Andy -- I'm not 100% sure why this is necessary.
    numpy.seterr(all='ignore')
    warnings.simplefilter("ignore", RuntimeWarning)

    logger.info('CoverView v1.2.0 started running')


def create_gui_output_directory():
    _logger("Creating output directory structure for GUI")
    dir = options.output + '_gui'

    if os.path.exists(dir):
        shutil.rmtree(dir)

    os.makedirs(dir)
    os.makedirs(dir + '/data')
    cvdir = os.path.dirname(os.path.realpath(__file__))
    shutil.copy(cvdir + '/gui.html', dir + '/' + options.output + '_coverview.html')
    shutil.copytree(cvdir + '/lib', dir + '/lib')


if __name__ == "__main__":

    configure_logging()
    options, config = get_input_options()

    _logger.info("Running CoverView {} with options".format(__version))
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

        _logger.info('CoverView {} succesfully finished'.format(__version))
    else:
        if config['outputs']['gui']:
            create_gui_output_directory()

        target_names, unique_target_ids = get_names_of_target_regions(options.bedfile)
        number_of_targets = len(target_names)
        sample_name = options.input.replace(".bam", "")

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
            output.finalizeJSONOutput(
                options,
                chromosome_coverage_metrics,
                config,
                number_of_targets,
                num_failed_targets,
                unique_target_ids,
                ids_of_failed_targets
            )

        _logger.info("CoverView {} succesfully finished".format(__version))

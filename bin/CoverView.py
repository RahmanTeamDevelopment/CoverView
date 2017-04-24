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

from coverview import transcript
from coverview import output
from coverview.calculators import calculateChromData, get_profiles


logger = logging.getLogger("coverview")


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


class SingleJob(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.samfile = pysam.Samfile(options.input, "rb")
        self.entsdb = None
        self.out_targets = None
        self.out_poor = None
        self.out_json = None
        self.out_profiles = None
        self.num_reads_on_target = {}
        self.ids_of_failed_targets = set()

        if config['outputs']['gui']:
            self.reffile = pysam.Fastafile(config["reference"])

        if not config['transcript_db'] is None:
            self.enstdb = pysam.Tabixfile(config['transcript_db'])

        if config['outputs']['regions']:
            self.out_targets = open(options.output + '_regions.txt', 'w')

        if config['outputs']['profiles']:
            self.out_profiles = open(options.output + '_profiles.txt', 'w')
            if not config['transcript_db'] is None:
                self.out_poor = open(options.output + '_poor.txt', 'w')

        if config['outputs']['gui']:
            self.out_json = open(options.output + '_gui/data/results.js', 'w')

        self.first = True

    def targetPASS(self, target):
        summary = target['Summary']
        for key, value in self.config['pass'].iteritems():
            [metrics, minmax] = key.split('_')
            if minmax == 'MIN' and summary[str(metrics)] < float(value):
                return False
            if minmax == 'MAX' and summary[str(metrics)] > float(value):
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
        if self.config['outputs']['regions']:
            self.out_targets.close()

        if self.config['outputs']['profiles']:
            self.out_profiles.close()
            if not config['transcript_db'] is None:
                self.out_poor.close()

        if self.config['outputs']['gui']:
            self.out_json.write(']')
            self.out_json.close()

    def run(self):
        logger.info("Coverage metrics will be generated in a single process")
        logger.info("Writing output headers")

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
            for target in get_profiles(self.samfile, cluster, self.config):

                if target is None:
                    continue

                chrom = target['Chrom']
                begin = target['Start']
                end = target['End']
                self.num_reads_on_target[target['Name']] = target['Profiles']['RC']

                if config['outputs']['regions']:
                    target['Summary'] = self.compute_summaries_of_region_coverage(target['Profiles'])

                    if config['pass'] is not None:
                        target['PASS'] = self.targetPASS(target)
                    else:
                        target['PASS'] = True

                    if self.config['outputs']['gui']:
                        target['Ref'] = self.getReferenceSequence(chrom, begin, end)
                    else:
                        target['Ref'] = ''
                else:
                    target['PASS'] = True

                if not target['PASS']:
                    if '_' in target['Name']:
                        ids = target['Name'][:target['Name'].find('_')]
                    else:
                        ids = target['Name']

                    self.ids_of_failed_targets.append(ids)
                self.output(target)

        self.close_output_files()
        logger.info("Finished computing coverage metrics in all regions")
        logger.info("Data was processed in {} clusters".format(num_clusters))

    def getReferenceSequence(self, chrom, start, end):
        goodchrom = chrom
        if not goodchrom in self.reffile.references:
            goodchrom = 'chr' + chrom
            if not goodchrom in self.reffile.references:
                goodchrom = chrom[3:]
                if not goodchrom in self.reffile.references:
                    return None

        if start < 0:
            start = 1

        last = self.reffile.getReferenceLength(goodchrom)

        if end > last:
            end = last

        seq = self.reffile.fetch(goodchrom, start - 1, end)
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

        # Calculate transcript coordinates for the _regions output file
        if self.config['transcript']['regions'] and not self.config['transcript_db'] is None:

            transcoords_start = transcript.getTranscriptCoordinates(self.enstdb, target['Chrom'], target['Start'])
            transcoords_end = transcript.getTranscriptCoordinates(self.enstdb, target['Chrom'], target['End'])

            v = []

            for key, value in transcoords_start.iteritems():
                v.append(key.geneSymbol + ':' + key.ENST + ':' + value)

            transcoordstr_start = ','.join(v)
            v = []

            for key, value in transcoords_end.iteritems():
                v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
            transcoordstr_end = ','.join(v)

        record = [target['Name'], target['Chrom'], str(target['Start']), str(target['End'])]

        if self.config['transcript']['regions'] and not self.config['transcript_db'] is None:
            if not transcoordstr_start == '':
                record.extend([transcoordstr_start, transcoordstr_end])
            else:
                record.extend(['.', '.'])

        if not self.config['pass'] is None:
            if target['PASS']: record.append('PASS')
            else: record.append('FAIL')

        if str(summary['MAXFLMQ']) == 'nan': summary['MAXFLMQ'] = '.'
        if str(summary['MAXFLBQ']) == 'nan': summary['MAXFLBQ'] = '.'

        record.extend([str(target['Profiles']['RC']), str(summary['MEDCOV']), str(summary['MINCOV']), str(summary['MEDQCOV']),
                       str(summary['MINQCOV']), str(summary['MAXFLMQ']), str(summary['MAXFLBQ'])])

        if self.config['direction']:
            if str(summary['MAXFLMQ_f']) == 'nan': summary['MAXFLMQ_f'] = '.'
            if str(summary['MAXFLBQ_f']) == 'nan': summary['MAXFLBQ_f'] = '.'
            if str(summary['MAXFLMQ_r']) == 'nan': summary['MAXFLMQ_r'] = '.'
            if str(summary['MAXFLBQ_r']) == 'nan': summary['MAXFLBQ_r'] = '.'

            record.extend([str(summary['MEDCOV_f']), str(summary['MINCOV_f']), str(summary['MEDQCOV_f']),
                           str(summary['MINQCOV_f']), str(summary['MAXFLMQ_f']), str(summary['MAXFLBQ_f'])])

            record.extend([str(summary['MEDCOV_r']), str(summary['MINCOV_r']), str(summary['MEDQCOV_r']),
                           str(summary['MINQCOV_r']), str(summary['MAXFLMQ_r']), str(summary['MAXFLBQ_r'])])

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


if __name__ == "__main__":

    logger = logging.getLogger("coverview")

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(pathname)s - Line %(lineno)s - %(message)s")

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)

    logger.addHandler(stream_handler)
    logger.setLevel(logging.DEBUG)

    numpy.seterr(all='ignore')
    warnings.simplefilter("ignore", RuntimeWarning)

    logger.info('CoverView v1.2.0 started running')

    options, config = get_input_options()

    logger.info("Configuration options")
    logger.info(options)
    logger.info(config)

    if options.bedfile is None:
        logger.info("No input BED file specified. Computing minimal coverage information")
        samfile = pysam.Samfile(options.input, "rb")
        chromdata = calculateChromdata_minimal(samfile, options)
        output.output_summary_minimal(options, chromdata)
        logger.info('CoverView v1.2.0 succesfully finished')
    else:
        if config['outputs']['gui']:
            dir = options.output + '_gui'

            if os.path.exists(dir):
                shutil.rmtree(dir)

            os.makedirs(dir)
            os.makedirs(dir+'/data')
            cvdir = os.path.dirname(os.path.realpath(__file__))
            shutil.copy(cvdir+'/gui.html', dir+'/'+options.output+'_coverview.html')
            shutil.copytree(cvdir+'/lib', dir+'/lib')

        names, uniqueIDs = makeNames(options.bedfile)
        number_of_targets = len(names)
        samplename = options.input[:options.input.rfind('.bam')]
        logger.info("There are {} target regions".format(number_of_targets))

        process = SingleJob(options, config)
        process.run()

        ontarget = process.num_reads_on_target

        ids_of_failed_targets = process.ids_of_failed_targets
        num_failed_targets = len(ids_of_failed_targets)

        logger.info("{} regions failed the coverage thresholds".format(num_failed_targets))

        samfile = pysam.Samfile(options.input, "rb")
        chromdata = calculateChromData(samfile, ontarget)
        output.output_summary(options, chromdata)

        if config['outputs']['gui']:
            output.finalizeJSONOutput(
                options,
                chromdata,
                config,
                number_of_targets,
                num_failed_targets,
                uniqueIDs,
                ids_of_failed_targets
            )

        logger.info("CoverView v1.2.0 succesfully finished")

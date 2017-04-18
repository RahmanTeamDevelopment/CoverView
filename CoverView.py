#!/env/bin/python

from __future__ import division

from optparse import OptionParser

import coverage
import datetime
import json
import numpy
import os
import pysam
import shutil
import sys
import transcript
import warnings

from coverage.output import *


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


class SingleJob(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.reads = dict()
        self.samfile = pysam.Samfile(options.input, "rb")

        if config['outputs']['gui']:
            self.reffile = pysam.Fastafile(config["reference"])

        if not config['transcript_db'] is None:
            self.enstdb = pysam.Tabixfile(config['transcript_db'])
        else:
            self.enstdb = None

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

    def get_profiles(self, region):
        """
        Calculate and return coverage metrics for a specified region. Metrics include total
        coverage, coverage above the required base-quality and mapping quality threshold, fractions of
        low base qualities at a given position, fractions of low mapping quality reads covering a given
        position, median mapping qualities and median base qualities.

        This is by far the most computationally expensive part of CoverView. > 90% of the run-time is
        currently spent in this function.
        """
        return coverage.get_profiles(self.samfile, region, self.config)

    def run(self):
        self.run_process()

    def run_process(self):
        output_target_file_header(
            self.config,
            self.out_poor,
            self.out_json,
            self.out_targets,
            self.out_profiles
        )

        print ''
        sys.stdout.write('\rRunning analysis ... 0.0%')
        sys.stdout.flush()

        fails = []
        numOfFails = 0

        for index,line in enumerate(open(self.options.bedfile)):
            target = dict()
            line = line.rstrip()

            if line == '':
                continue

            if line.startswith('#'):
                continue

            row = line.split('\t')

            if len(row) < 4:
                print "Error: incorrect BED file!"
                quit()

            chrom = row[0]
            begin = int(row[1])
            end = int(row[2])
            key = row[3]

            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            region = chrom + ':' + row[1] + '-' + row[2]

            target['Name'] = key
            target['Chrom'] = chrom
            target['Start'] = begin + 1
            target['End'] = end

            summary = dict()

            if self.config['outputs']['profiles'] or self.config['outputs']['regions']:
                profiles = dict()

                if self.config['direction']:
                    count, count_f, count_r,COV, QCOV, MEDBQ, FLBQ, MEDMQ, FLMQ, COV_f, QCOV_f, MEDBQ_f, FLBQ_f, MEDMQ_f, FLMQ_f, COV_r, QCOV_r, MEDBQ_r, FLBQ_r, MEDMQ_r, FLMQ_r = self.get_profiles(region)
                    summary['RC'] = count
                    summary['RC_f'] = count_f
                    summary['RC_r'] = count_r
                else:
                    count, COV, QCOV, MEDBQ, FLBQ, MEDMQ, FLMQ = self.get_profiles(region)
                    summary['RC'] = count

                profiles['COV'] = COV
                profiles['QCOV'] = QCOV
                profiles['MEDBQ'] = MEDBQ
                profiles['FLBQ'] = FLBQ
                profiles['MEDMQ'] = MEDMQ
                profiles['FLMQ'] = FLMQ

                if self.config['direction']:
                    profiles['COV_f'] = COV_f
                    profiles['QCOV_f'] = QCOV_f
                    profiles['MEDBQ_f'] = MEDBQ_f
                    profiles['FLBQ_f'] = FLBQ_f
                    profiles['MEDMQ_f'] = MEDMQ_f
                    profiles['FLMQ_f'] = FLMQ_f
                    profiles['COV_r'] = COV_r
                    profiles['QCOV_r'] = QCOV_r
                    profiles['MEDBQ_r'] = MEDBQ_r
                    profiles['FLBQ_r'] = FLBQ_r
                    profiles['MEDMQ_r'] = MEDMQ_r
                    profiles['FLMQ_r'] = FLMQ_r

                target['Profiles'] = profiles

            if self.config['outputs']['regions']:

                summary['MEDCOV'] = numpy.median(COV)
                summary['MEDQCOV'] = numpy.median(QCOV)
                summary['MINCOV'] = min_or_nan(COV)
                summary['MINQCOV'] = min_or_nan(QCOV)
                summary['MAXFLBQ'] = max_or_nan(FLBQ)
                summary['MAXFLMQ'] = max_or_nan(FLMQ)

                if self.config['direction']:

                    summary['MEDCOV_f'] = numpy.median(COV_f)
                    summary['MEDQCOV_r'] = numpy.median(QCOV_r)
                    summary['MEDQCOV_f'] = numpy.median(QCOV_f)
                    summary['MEDCOV_r'] = numpy.median(COV_r)
                    summary['MINCOV_f'] = min_or_nan(COV_f)
                    summary['MINQCOV_f'] = min_or_nan(QCOV_f)
                    summary['MAXFLBQ_f'] = max_or_nan(FLBQ_f)
                    summary['MAXFLMQ_f'] = max_or_nan(FLMQ_f)
                    summary['MINCOV_r'] = min_or_nan(COV_r)
                    summary['MINQCOV_r'] = min_or_nan(QCOV_r)
                    summary['MAXFLBQ_r'] = max_or_nan(FLBQ_r)
                    summary['MAXFLMQ_r'] = max_or_nan(FLMQ_r)

                # Add summary metrics to target dict
                target['Summary'] = summary

                # Calculate if target is PASS or FAIL
                if not config['pass'] is None:
                    target['PASS'] = self.targetPASS(target)
                else:
                    target['PASS'] = True

                # Retrieve reference sequence if GUI output is created
                if self.config['outputs']['gui']:
                    target['Ref'] = self.getReferenceSequence(chrom, begin, end)
                else:
                    target['Ref'] = ''

            else:
                target['PASS'] = True

            self.output(target)

            if not target['PASS']:
                numOfFails += 1
                if '_' in target['Name']: ids = target['Name'][:target['Name'].find('_')]
                else: ids = target['Name']
                fails.append(ids)

            if index % 100 == 0:
                x = round(100 * index / numOfTargets, 1)
                x = min(x, 100.0)
                sys.stdout.write('\rRunning analysis ... ' + str(x) + '%')
                sys.stdout.flush()

        if self.config['outputs']['regions']:
            self.out_targets.close()

        if self.config['outputs']['profiles']:
            self.out_profiles.close()
            if not config['transcript_db'] is None:
                self.out_poor.close()

        if self.config['outputs']['gui']:
            self.out_json.write(']')
            self.out_json.close()

        with open(options.output + '_failedtargets' + '.txt', 'w') as failedtargetsfile:
            failedtargetsfile.write(str(fails) + '\n')

        with open(options.output + '_reads_on_target' + '.txt', 'w') as ontargetfile:
            for k, v in self.reads.iteritems():
                ontargetfile.write(k + ':' + str(len(v)) + '\n')

        sys.stdout.write('\rRunning analysis ... 100.0%')
        sys.stdout.flush()

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
                coverage.output.output_profiles(
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

        record.extend([str(summary['RC']), str(summary['MEDCOV']), str(summary['MINCOV']), str(summary['MEDQCOV']),
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


def defaultConfigs(config):
    if not 'duplicates' in config.keys(): config['duplicates'] = True
    if not 'outputs' in config.keys(): config['outputs'] = {'regions': True,
                                                            'profiles': True, 'gui': False}
    if not 'regions' in config['outputs'].keys(): config['outputs']['regions'] = True
    if not 'profiles' in config['outputs'].keys(): config['outputs']['profiles'] = True
    if not 'gui' in config['outputs'].keys(): config['outputs']['gui'] = False
    if not 'low_bq' in config.keys(): config['low_bq'] = 10
    if not 'low_mq' in config.keys(): config['low_mq'] = 20
    if not 'only_fail_profiles' in config.keys(): config['only_fail_profiles'] = False
    if not 'transcript_db' in config.keys(): config['transcript_db'] = None
    if not 'pass' in config.keys() or len(config['pass']) == 0: config['pass'] = None
    if not 'transcript' in config.keys(): config['transcript'] = {'regions': True, 'profiles': True}
    if not 'regions' in config['transcript'].keys(): config['transcript']['regions'] = True
    if not 'profiles' in config['transcript'].keys(): config['transcript']['profiles'] = True
    if not 'direction' in config.keys(): config['direction'] = False


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
    config = {}

    if options.config is not None:
        with open(options.config) as config_file:
            config = json.load(config_file)

    defaultConfigs(config)
    return options, config


if __name__ == "__main__":
    numpy.seterr(all='ignore')
    warnings.simplefilter("ignore", RuntimeWarning)

    print ""
    print "======================================================================================"
    print 'CoverView v1.2.0 started running: ', datetime.datetime.now()
    print ""

    options, config = get_input_options()

    if options.bedfile is None:
        printInfo_minimal(options)
        samfile = pysam.Samfile(options.input, "rb")
        chromdata = calculateChromdata_minimal(samfile)
        output_summary_minimal(options, chromdata)
        print ""
        print 'CoverView v1.2.0 succesfully finished: ', datetime.datetime.now()
        print "======================================================================================"
        print ""
        quit()

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
    numOfTargets = len(names)
    samplename = options.input[:options.input.rfind('.bam')]
    printInfo(options, config, numOfTargets)

    process = SingleJob(options, config)
    process.run()

    ontarget = dict()

    for line in open(options.output + '_reads_on_target' + '.txt'):
        [key, value] = line.split(':')
        if key in ontarget.keys():
            ontarget[key] += int(value.strip())
        else:
            ontarget[key] = int(value.strip())

    os.remove(options.output + '_reads_on_target' + '.txt')

    failedtargets = 0
    uniqueids = set()

    for line in open(options.output + '_failedtargets' + '.txt'):
        line = line.strip()
        fails = line[1:-1].split(',')
        failedtargets += len(fails)
        for x in fails:
            x = x.strip()
            uniqueids.add(x[1:-1])

    os.remove(options.output + '_failedtargets' + '.txt')

    print ' - Done. (' + str(failedtargets) + ' failed regions)'

    samfile = pysam.Samfile(options.input, "rb")
    chromdata = calculateChromdata(samfile, ontarget)
    output_summary(options, chromdata)

    if config['outputs']['gui']:
        finalizeJSONOutput(options)

    print ""
    print 'CoverView v1.2.0 succesfully finished: ', datetime.datetime.now()
    print "======================================================================================"
    print ""

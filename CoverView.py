#!/usr/bin/python

from __future__ import division
from optparse import OptionParser
import os
import sys
import datetime

# Checking if installation is complete
if not os.path.isfile(os.path.dirname(os.path.realpath(__file__))+"/pysamdir/pysam/Pileup.py"):
    print '\nCoverView installation not complete.'
    print 'Please run install.sh in the CoverView directory.\n'
    quit()

# Pysam import
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/pysamdir')
import pysam

import numpy
import multiprocessing
import json
import transcript
import warnings
import shutil


#########################################################################################################################################

class SingleJob(multiprocessing.Process):
    # Process constructor
    def __init__(self, threadidx, options, config, startline, endline, names):
        multiprocessing.Process.__init__(self)

        # Initializing variables
        self.threadidx = threadidx
        self.options = options
        self.config = config
        self.startline = startline
        self.endline = endline
        self.names = names

        # Initializing reads directory
        self.reads = dict()

        # Connecting to BAM file
        self.samfile = pysam.Samfile(options.input, "rb")

        # Connecting to reference genome
        if config['outputs']['gui']:
            self.reffile = pysam.Fastafile(config["reference"])

        # Connecting to transcript database file
        if not config['transcript_db'] is None:
            self.enstdb = pysam.Tabixfile(config['transcript_db'])
        else:
            self.enstdb = None

        # Initializing output files
        if int(options.threads) > 1:
            if config['outputs']['regions']:
                self.out_targets = open(options.output + '_regions_tmp_' + str(threadidx) + '.txt', 'w')
            if config['outputs']['profiles']:
                self.out_profiles = open(options.output + '_profiles_tmp_' + str(threadidx) + '.txt', 'w')
                if not config['transcript_db'] is None:
                    self.out_poor = open(options.output + '_poor_tmp_' + str(threadidx) + '.txt', 'w')
            if config['outputs']['gui']:
                self.out_json = open(options.output + '_results_tmp_' + str(threadidx) + '.json', 'w')
        else:
            if config['outputs']['regions']:
                self.out_targets = open(options.output + '_regions.txt', 'w')
            if config['outputs']['profiles']:
                self.out_profiles = open(options.output + '_profiles.txt', 'w')
                if not config['transcript_db'] is None:
                    self.out_poor = open(options.output + '_poor.txt', 'w')
            if config['outputs']['gui']:
                self.out_json = open(options.output + '_gui/data/results.js', 'w')

        self.first = True

    # Checking if target is fail or pass
    def targetPASS(self, target):
        summary = target['Summary']
        for key, value in self.config['pass'].iteritems():
            [metrics, minmax] = key.split('_')
            if minmax == 'MIN' and summary[str(metrics)] < float(value): return False
            if minmax == 'MAX' and summary[str(metrics)] > float(value): return False
        return True

    # Calculating COV, QCOV, MEDBQ, FLBQ, MEDMQ and FLMQ profiles for a region
    def getProfiles(self, region):

        # Base quality and mapping quality cutoff parameters
        bq_cutoff = float(self.config['low_bq'])
        mq_cutoff = float(self.config['low_mq'])

        # Splitting region to chrom:begin-end
        chrom = region[:region.find(':')]
        begin = int(region[region.find(':') + 1:region.find('-')]) - 1
        end = int(region[region.find('-') + 1:]) - 1

        # Initializing profiles
        ret_COV = [0] * (end - begin)
        ret_QCOV = [0] * (end - begin)
        ret_MEDBQ = [float('NaN')] * (end - begin)
        ret_FLBQ = [float('NaN')] * (end - begin)
        ret_MEDMQ = [float('NaN')] * (end - begin)
        ret_FLMQ = [float('NaN')] * (end - begin)

        if self.config['direction']:
            ret_COV_f = [0] * (end - begin)
            ret_QCOV_f = [0] * (end - begin)
            ret_MEDBQ_f = [float('NaN')] * (end - begin)
            ret_FLBQ_f = [float('NaN')] * (end - begin)
            ret_MEDMQ_f = [float('NaN')] * (end - begin)
            ret_FLMQ_f = [float('NaN')] * (end - begin)

            ret_COV_r = [0] * (end - begin)
            ret_QCOV_r = [0] * (end - begin)
            ret_MEDBQ_r = [float('NaN')] * (end - begin)
            ret_FLBQ_r = [float('NaN')] * (end - begin)
            ret_MEDMQ_r = [float('NaN')] * (end - begin)
            ret_FLMQ_r = [float('NaN')] * (end - begin)

        goodchrom = chrom
        chrprefix = self.samfile.references[0].startswith('chr')
        if chrprefix and not chrom.startswith('chr'): goodchrom = 'chr' + chrom
        if not chrprefix and chrom.startswith('chr'): goodchrom = chrom[3:]

        if not goodchrom in self.samfile.references: return None

        # Getting pileup for the given region
        if config['duplicates']: x = self.samfile.pileup(goodchrom, begin, end + 1, mask=0)
        else: x = self.samfile.pileup(goodchrom, begin, end + 1)

        i = 0
        # Iterating through columns of the pileup
        for pileupcolumn in x:

            if pileupcolumn.pos > begin and pileupcolumn.pos <= end:

                # Getting coverage
                ret_COV[i] = pileupcolumn.n

                if self.config['direction']:
                    cov_f = 0
                    cov_r = 0
                    bqs_f = []
                    bqs_r = []
                    mqs_f = []
                    mqs_r = []
                    qcov_f = 0
                    qcov_r = 0

                bqs = []
                mqs = []
                qcov = 0
                # Iterating through reads mapping to the pileup column
                for pileupread in pileupcolumn.pileups:

                    if self.config['direction']:
                        if pileupread.alignment.is_reverse: cov_r += 1
                        else: cov_f += 1

                    # Checking if it is a deletion
                    if pileupread.is_del:
                        if pileupread.alignment.mapq >= mq_cutoff:
                            qcov += 1
                            if self.config['direction']:
                                if pileupread.alignment.is_reverse:
                                    qcov_r += 1
                                else:
                                    qcov_f += 1
                        continue


                    # Getting base quality
                    c = pileupread.alignment.qual[pileupread.qpos]
                    bq = ord(c) - 33
                    bqs.append(bq)

                    if self.config['direction']:
                        if pileupread.alignment.is_reverse: bqs_r.append(bq)
                        else: bqs_f.append(bq)

                    # Getting mapping quality
                    mqs.append(pileupread.alignment.mapq)

                    if self.config['direction']:
                        if pileupread.alignment.is_reverse: mqs_r.append(pileupread.alignment.mapq)
                        else: mqs_f.append(pileupread.alignment.mapq)

                    # Calculating QCOV
                    if pileupread.alignment.mapq >= mq_cutoff and bq >= bq_cutoff:
                        qcov += 1
                        if self.config['direction']:
                            if pileupread.alignment.is_reverse: qcov_r += 1
                            else: qcov_f += 1

                # Summary stats of QCOV, MEDBQ, FLBQ, MEDMQ and FLMQ
                ret_QCOV[i] = qcov
                ret_MEDBQ[i] = numpy.median(bqs)
                ret_MEDMQ[i] = numpy.median(mqs)
                if len(bqs) > 0: ret_FLBQ[i] = round(len([x for x in bqs if x < bq_cutoff]) / len(bqs), 3)
                if len(mqs) > 0: ret_FLMQ[i] = round(len([x for x in mqs if x < mq_cutoff]) / len(mqs), 3)

                # Directionality
                if self.config['direction']:
                    ret_COV_f[i] = cov_f
                    ret_COV_r[i] = cov_r
                    ret_QCOV_f[i] = qcov_f
                    ret_QCOV_r[i] = qcov_r
                    ret_MEDBQ_f[i] = numpy.median(bqs_f)
                    ret_MEDBQ_r[i] = numpy.median(bqs_r)
                    ret_MEDMQ_f[i] = numpy.median(mqs_f)
                    ret_MEDMQ_r[i] = numpy.median(mqs_r)

                    if len(bqs_f) > 0: ret_FLBQ_f[i] = round(len([x for x in bqs_f if x < bq_cutoff]) / len(bqs_f), 3)
                    if len(bqs_r) > 0: ret_FLBQ_r[i] = round(len([x for x in bqs_r if x < bq_cutoff]) / len(bqs_r), 3)
                    if len(mqs_f) > 0: ret_FLMQ_f[i] = round(len([x for x in mqs_f if x < mq_cutoff]) / len(mqs_f), 3)
                    if len(mqs_r) > 0: ret_FLMQ_r[i] = round(len([x for x in mqs_r if x < mq_cutoff]) / len(mqs_r), 3)

                i += 1

        if self.config['direction']: return ret_COV, ret_QCOV, ret_MEDBQ, ret_FLBQ, ret_MEDMQ, ret_FLMQ, ret_COV_f, ret_QCOV_f, ret_MEDBQ_f, ret_FLBQ_f, ret_MEDMQ_f, ret_FLMQ_f, ret_COV_r, ret_QCOV_r, ret_MEDBQ_r, ret_FLBQ_r, ret_MEDMQ_r, ret_FLMQ_r
        else: return ret_COV, ret_QCOV, ret_MEDBQ, ret_FLBQ, ret_MEDMQ, ret_FLMQ

    # Calculating read counts for a region
    def readcountsForRegion(self, region):

        # Splitting region to chrom:begin-end
        chrom = region[:region.find(':')]
        begin = int(region[region.find(':') + 1:region.find('-')])
        end = int(region[region.find('-') + 1:])

        if not chrom in self.reads.keys(): self.reads[chrom] = set()

        goodchrom = chrom
        chrprefix = self.samfile.references[0].startswith('chr')
        if chrprefix and not chrom.startswith('chr'): goodchrom = 'chr' + chrom
        if not chrprefix and chrom.startswith('chr'): goodchrom = chrom[3:]

        if not goodchrom in self.samfile.references: return None, None, None

        count = 0
        count_f = 0
        count_r = 0
        for x in self.samfile.fetch(goodchrom, begin, end):
            if self.config['duplicates']:
                count += 1
                if x.is_reverse: count_r += 1
                else: count_f += 1
            else:
                if not x.is_duplicate:
                    count += 1
                    if x.is_reverse: count_r += 1
                    else: count_f += 1
            self.reads[chrom].add((x.qname, x.pos))

        return count, count_f, count_r

    # Calculating read counts for a region
    def readcountsForRegionOLD(self, region):

        # Splitting region to chrom:begin-end
        chrom = region[:region.find(':')]
        begin = int(region[region.find(':') + 1:region.find('-')])
        end = int(region[region.find('-') + 1:])

        if not chrom in self.reads.keys(): self.reads[chrom] = set()

        goodchrom = chrom
        chrprefix = self.samfile.references[0].startswith('chr')
        if chrprefix and not chrom.startswith('chr'): goodchrom = 'chr' + chrom
        if not chrprefix and chrom.startswith('chr'): goodchrom = chrom[3:]

        if not goodchrom in self.samfile.references: return None

        count = 0
        for x in self.samfile.fetch(goodchrom, begin, end):
            if self.config['duplicates']:
                count += 1
            else:
                if not x.is_duplicate: count += 1
            self.reads[chrom].add((x.qname, x.pos))

        return count

    # Running process
    def run(self):

        # If there is only one thread
        if int(options.threads) == 1:

            # Outputting _target file header
            if self.config['outputs']['regions']:
                targetheader = ['Region', 'Chromosome', 'Start_position', 'End_position']
                if config['transcript']['regions'] and not config['transcript_db'] is None: targetheader.extend(
                    ['Start_transcript', 'End_transcript'])
                if not config['pass'] is None: targetheader.append('Pass_or_fail')
                targetheader.extend(['RC', 'MEDCOV', 'MINCOV', 'MEDQCOV', 'MINQCOV', 'MAXFLMQ', 'MAXFLBQ'])

                if self.config['direction']:
                    targetheader.extend(['MEDCOV+', 'MINCOV+', 'MEDQCOV+', 'MINQCOV+', 'MAXFLMQ+', 'MAXFLBQ+'])
                    targetheader.extend(['MEDCOV-', 'MINCOV-', 'MEDQCOV-', 'MINQCOV-', 'MAXFLMQ-', 'MAXFLBQ-'])

                self.out_targets.write('#' + '\t'.join(targetheader) + '\n')

            # Outputting _profiles file header
            if self.config['outputs']['profiles']:
                profheader = ['Chromosome', 'Position']
                if config['transcript']['profiles'] and not config['transcript_db'] is None: profheader.append(
                    'Transcript_coordinate')
                profheader.extend(['COV', 'QCOV', 'MEDBQ', 'FLBQ', 'MEDMQ', 'FLMQ'])

                if self.config['direction']:
                    profheader.extend(['COV+', 'QCOV+', 'MEDBQ+', 'FLBQ+', 'MEDMQ+', 'FLMQ+'])
                    profheader.extend(['COV-', 'QCOV-', 'MEDBQ-', 'FLBQ-', 'MEDMQ-', 'FLMQ-'])

                self.out_profiles.write('#' + '\t'.join(profheader) + '\n')

            # Outputting _poor file header
            if not config['transcript_db'] is None and config['outputs']['profiles']:
                poorheader = ['Region', 'Chromosome', 'Start_position', 'End_position', 'Start_transcript',
                              'End_transcript']
                self.out_poor.write('#' + '\t'.join(poorheader) + '\n')

            # Initializing gui data json file
            if self.config['outputs']['gui']:
                self.out_json.write('function readData() {\n')
                self.out_json.write('\tdata={\"targets\":[')

        # Initializing progress info
        if self.threadidx == 1:
            print ''
            sys.stdout.write('\rRunning analysis ... 0.0%')
            sys.stdout.flush()

        # Iteriate through the bed file
        fails = []
        numOfFails = 0
        counter = 0
        for line in open(self.options.bedfile):
            counter += 1

            # Check which part of the bed file should be process by this thread
            if counter < int(self.startline): continue
            if not self.endline == '':
                if counter > int(self.endline): break

            # Retrieve target as a dictionary
            target = dict()

            line = line.rstrip()
            if line == '': continue
            if line.startswith('#'): continue

            row = line.split('\t')
            if len(row) < 4:
                print "Error: incorrect BED file!"
                quit()

            key = self.names[counter - 1]

            chrom = row[0]
            begin = int(row[1])
            end = int(row[2])
            if not chrom.startswith('chr'): chrom = 'chr' + chrom
            region = chrom + ':' + row[1] + '-' + row[2]

            target['Name'] = key
            target['Chrom'] = chrom
            target['Start'] = begin + 1
            target['End'] = end

            # Calculate profiles if either _profiles or _regions files are to be outputted
            if self.config['outputs']['profiles'] or self.config['outputs']['regions']:
                profiles = dict()

                # Calculate directional and/or non-directional profiles
                if self.config['direction']: COV, QCOV, MEDBQ, FLBQ, MEDMQ, FLMQ, COV_f, QCOV_f, MEDBQ_f, FLBQ_f, MEDMQ_f, FLMQ_f, COV_r, QCOV_r, MEDBQ_r, FLBQ_r, MEDMQ_r, FLMQ_r = self.getProfiles(region)
                else: COV, QCOV, MEDBQ, FLBQ, MEDMQ, FLMQ = self.getProfiles(region)

                # Make profiles dictionary
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

                # Add profiles to target dict
                target['Profiles'] = profiles


            # Calculate summary metrics for the target

            # Read counts
            summary = dict()
            summary['RC'], count_f, count_r = self.readcountsForRegion(region)
            if self.config['direction']:
                summary['RC_f'] = count_f
                summary['RC_r'] = count_r

            # Calculate additional summary metrics if _regions file is to be outputted
            if self.config['outputs']['regions']:

                # Caclulate MEDCOV, MINCOV, MEDQCOV, MINQCOV, MAXFLBQ, MAXFLMQ
                summary['MEDCOV'] = numpy.median(COV)
                if len(COV) > 0: summary['MINCOV'] = min(COV)
                else: summary['MINCOV'] = float('NaN')
                summary['MEDQCOV'] = numpy.median(QCOV)
                if len(QCOV) > 0: summary['MINQCOV'] = min(QCOV)
                else: summary['MINQCOV'] = float('NaN')
                if len(FLBQ) > 0: summary['MAXFLBQ'] = max(FLBQ)
                else: summary['MAXFLBQ'] = float('NaN')
                if len(FLMQ) > 0: summary['MAXFLMQ'] = max(FLMQ)
                else: summary['MAXFLMQ'] = float('NaN')

                # Caclulate directional MEDCOV, MINCOV, MEDQCOV, MINQCOV, MAXFLBQ, MAXFLMQ
                if self.config['direction']:

                    summary['MEDCOV_f'] = numpy.median(COV_f)
                    if len(COV_f) > 0: summary['MINCOV_f'] = min(COV_f)
                    else: summary['MINCOV_f'] = float('NaN')
                    summary['MEDQCOV_f'] = numpy.median(QCOV_f)
                    if len(QCOV_f) > 0: summary['MINQCOV_f'] = min(QCOV_f)
                    else: summary['MINQCOV_f'] = float('NaN')
                    if len(FLBQ_f) > 0: summary['MAXFLBQ_f'] = max(FLBQ_f)
                    else: summary['MAXFLBQ_f'] = float('NaN')
                    if len(FLMQ_f) > 0: summary['MAXFLMQ_f'] = max(FLMQ_f)
                    else: summary['MAXFLMQ_f'] = float('NaN')

                    summary['MEDCOV_r'] = numpy.median(COV_r)
                    if len(COV_r) > 0: summary['MINCOV_r'] = min(COV_r)
                    else: summary['MINCOV_r'] = float('NaN')
                    summary['MEDQCOV_r'] = numpy.median(QCOV_r)
                    if len(QCOV_r) > 0: summary['MINQCOV_r'] = min(QCOV_r)
                    else: summary['MINQCOV_r'] = float('NaN')
                    if len(FLBQ_r) > 0: summary['MAXFLBQ_r'] = max(FLBQ_r)
                    else: summary['MAXFLBQ_r'] = float('NaN')
                    if len(FLMQ_r) > 0: summary['MAXFLMQ_r'] = max(FLMQ_r)
                    else: summary['MAXFLMQ_r'] = float('NaN')

                # Add summary metrics to target dict
                target['Summary'] = summary

                # Calculate if target is PASS or FAIL
                if not config['pass'] is None: target['PASS'] = self.targetPASS(target)
                else: target['PASS'] = True

                # Retrieve reference sequence if GUI output is created
                if self.config['outputs']['gui']: target['Ref'] = self.getReferenceSequence(chrom, begin, end)
                else: target['Ref'] = ''

            else:
                target['PASS'] = True

            # Output target data
            self.output(target)

            # Count failed targets
            if not target['PASS']:
                numOfFails += 1
                if '_' in target['Name']: ids = target['Name'][:target['Name'].find('_')]
                else: ids = target['Name']
                fails.append(ids)

            # Update progress info
            if self.threadidx == 1:
                if counter % 100 == 0:
                    x = round(100 * counter / (numOfTargets / int(self.options.threads)), 1)
                    x = min(x, 100.0)
                    sys.stdout.write('\rRunning analysis ... ' + str(x) + '%')
                    sys.stdout.flush()

        # Close output files
        if self.config['outputs']['regions']: self.out_targets.close()
        if self.config['outputs']['profiles']:
            self.out_profiles.close()
            if not config['transcript_db'] is None:
                self.out_poor.close()

        if self.config['outputs']['gui']:
            if int(options.threads) == 1: self.out_json.write(']')
            self.out_json.close()

        # Write failed targets and reads on target to tmp files
        with open(options.output + '_failedtargets_' + str(self.threadidx) + '.txt', 'w') as failedtargetsfile:
            failedtargetsfile.write(str(fails) + '\n')

        with open(options.output + '_reads_on_target_' + str(self.threadidx) + '.txt', 'w') as ontargetfile:
            for k, v in self.reads.iteritems():
                ontargetfile.write(k + ':' + str(len(v)) + '\n')

        # Finalize progress info
        if self.threadidx == 1:
            sys.stdout.write('\rRunning analysis ... 100.0%')
            sys.stdout.flush()

    # Return reference sequence
    def getReferenceSequence(self, chrom, start, end):
        goodchrom = chrom
        if not goodchrom in self.reffile.references:
            goodchrom = 'chr' + chrom
            if not goodchrom in self.reffile.references:
                goodchrom = chrom[3:]
                if not goodchrom in self.reffile.references:
                    return None

        if start < 0: start = 1
        last = self.reffile.getReferenceLength(goodchrom)
        if end > last: end = last
        seq = self.reffile.fetch(goodchrom, start - 1, end)
        return seq.upper()

    # Outputting target
    def output(self, target):
        # _targets.txt file:
        if self.config['outputs']['regions']: self.output_target(target)

        # _profiles.txt and _poor.txt files:
        if self.config['outputs']['profiles']:
            if not self.config['only_fail_profiles']:
                self.output_profiles(target)
            else:
                if not self.config['pass'] is None:
                    if not target['PASS']:
                        self.output_profiles(target)
                else:
                    self.output_profiles(target)

        # Outputting json:
        if self.config['outputs']['gui']: self.output_json(target)
        if self.first: self.first = False

    # Writing _target output file
    def output_json(self, target):
        if not self.first: self.out_json.write(',')
        self.out_json.write(json.dumps(target, separators=(',', ':')))

    # Writing _target output file
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

        # Create output record
        record = [target['Name'], target['Chrom'], str(target['Start']), str(target['End'])]

        if self.config['transcript']['regions'] and not self.config['transcript_db'] is None:
            if not transcoordstr_start == '': record.extend([transcoordstr_start, transcoordstr_end])
            else: record.extend(['.', '.'])

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

        # Write record to output file
        self.out_targets.write('\t'.join(record) + '\n')

    # Writing _profiles output file
    def output_profiles(self, target):
        profiles = target['Profiles']

        # Write target name to output file
        self.out_profiles.write('\n')
        self.out_profiles.write('[' + target['Name'] + ']\n')

        if not self.config['transcript_db'] is None:
            window_qcov = {"targetname": None, "chrom": None, "start": None, "transcriptstart": None}

        # Iterate through all bases of the targeted region
        for i in range(len(profiles['COV'])):
            transcoordstr = ''

            # Calculate transcript coordinates if required
            if self.config['transcript']['profiles'] and not self.config['transcript_db'] is None:
                transcoords = transcript.getTranscriptCoordinates(self.enstdb, target['Chrom'], target['Start'] + i)
                v = []
                for key, value in transcoords.iteritems():
                    v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
                transcoordstr = ','.join(v)

            # Create output record
            record = [target['Chrom'], str(target['Start'] + i)]
            if self.config['transcript']['profiles'] and not self.config['transcript_db'] is None:
                record.append(transcoordstr)

            record.extend(
                [str(profiles['COV'][i]), str(profiles['QCOV'][i]), str(profiles['MEDBQ'][i]), str(profiles['FLBQ'][i]),
                 str(profiles['MEDMQ'][i]), str(profiles['FLMQ'][i])])

            if self.config['direction']:
                record.extend([str(profiles['COV_f'][i]), str(profiles['QCOV_f'][i]), str(profiles['MEDBQ_f'][i]),
                               str(profiles['FLBQ_f'][i]), str(profiles['MEDMQ_f'][i]), str(profiles['FLMQ_f'][i])])
                record.extend([str(profiles['COV_r'][i]), str(profiles['QCOV_r'][i]), str(profiles['MEDBQ_r'][i]),
                               str(profiles['FLBQ_r'][i]), str(profiles['MEDMQ_r'][i]), str(profiles['FLMQ_r'][i])])

            if record[-1] == 'nan': record[-1] = '.'
            if record[-2] == 'nan': record[-2] = '.'
            if record[-3] == 'nan': record[-3] = '.'
            if record[-4] == 'nan': record[-4] = '.'

            # Write record to output file
            self.out_profiles.write('\t'.join(record) + '\n')

            # Calculate poor quality windows and transcript coordinates for _poor output file
            if self.config['transcript']['poor'] and not self.config['transcript_db'] is None:
                if profiles['QCOV'][i] < 15:
                    if window_qcov['transcriptstart'] is None:
                        window_qcov['targetname'] = target['Name']
                        window_qcov['chrom'] = target['Chrom']
                        window_qcov['start'] = target['Start'] + i
                        if transcoordstr == '':
                            transcoords = transcript.getTranscriptCoordinates(self.enstdb, target['Chrom'],
                                                                              target['Start'] + i)
                            v = []
                            for key, value in transcoords.iteritems():
                                v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
                            transcoordstr = ','.join(v)
                        window_qcov['transcriptstart'] = transcoordstr
                else:
                    if not window_qcov['transcriptstart'] is None:

                        transcoords = transcript.getTranscriptCoordinates(self.enstdb, target['Chrom'],
                                                                          target['Start'] + i - 1)
                        v = []
                        for key, value in transcoords.iteritems():
                            v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
                        transcoordstr_end = ','.join(v)

                        record = [window_qcov['targetname'], window_qcov['chrom'], str(window_qcov['start']),
                                  str(target['Start'] + i - 1), window_qcov['transcriptstart'], transcoordstr_end]
                        if record[4] == '': record[4] = 'None'
                        if record[5] == '': record[5] = 'None'

                        if not (record[4] == 'None' and record[5] == 'None'):
                            self.out_poor.write('\t'.join(record) + '\n')

                        window_qcov = {"targetname": None, "chrom": None, "start": None, "transcriptstart": None}

        # Finalize calculation of poor quality windows and transcript coordinates for _poor output file
        if self.config['transcript']['poor'] and not self.config['transcript_db'] is None:
            if not window_qcov['transcriptstart'] is None:

                if transcoordstr == '':
                    transcoords = transcript.getTranscriptCoordinates(self.enstdb, target['Chrom'], target['Start'] + i)
                    v = []
                    for key, value in transcoords.iteritems():
                        v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
                    transcoordstr = ','.join(v)

                record = [window_qcov['targetname'], window_qcov['chrom'], str(window_qcov['start']),
                          str(target['Start'] + i), window_qcov['transcriptstart'], transcoordstr]
                if record[4] == '': record[4] = 'None'
                if record[5] == '': record[5] = 'None'

                if not (record[4] == 'None' and record[5] == 'None'):
                    self.out_poor.write('\t'.join(record) + '\n')

               # window_qcov = {"targetname": None, "chrom": None, "start": None, "transcriptstart": None}


#########################################################################################################################################

# Printing out info about parameters
def printInfo(options, config, numOfTargets):
    targetstxt = ' (' + str(numOfTargets) + ' regions)'
    print 'Input, output and settings:'
    print "--------------------------------------------------------------------------------------"

    if not options.config is None:
        print "Configuration file:     " + options.config
    else:
        print "Configuration file:     using default settings"

    print "Input file name:        " + options.input
    print "BED file name:          " + options.bedfile + targetstxt
    print ''
    if not config['transcript_db'] is None and (config['outputs']['regions'] or config['outputs']['profiles']) and (
            config['transcript']['regions'] or config['transcript']['profiles'] or config['transcript']['poor']):
        print "Transcript db file:     " + config['transcript_db']
        print ''

    formats = 'summary'
    if config['outputs']['regions']: formats += ', regions'
    if config['outputs']['profiles']:
        if config['only_fail_profiles']:
            formats += ', profiles (failed regions)'
        else:
            formats += ', profiles (all regions)'
    if not config['transcript_db'] is None and config['outputs']['profiles'] and config['transcript']['poor']:
        formats += ', poor'
    if config['outputs']['gui']: formats += ', GUI'
    print "Output formats:         " + formats
    print "Output files prefix:    " + options.output
    print ''
    if config['duplicates']:
        print "Duplicate reads:        Included"
    else:
        print "Duplicate reads:        Excluded"
    if config['direction']:
        print  "Directional metrics:    Outputted"
    else:
        print  "Directional metrics:    Not outputted"
    if config['outputs']['regions'] or config['outputs']['profiles']:
        print "Mapping quality cutoff: " + str(config['low_mq'])
        print "Base quality cutoff:    " + str(config['low_bq'])
    if not config['pass'] is None and (config['outputs']['regions'] or config['outputs']['profiles']):
        print ''
        params = []
        for k, v in config['pass'].iteritems():
            params.append(str(k) + '=' + str(v))
        print "Region pass parameters: " + '; '.join(params)
    if int(options.threads) > 1:
        print ''
        print 'Multithreading:         ' + str(options.threads) + ' processes'
    print "--------------------------------------------------------------------------------------"

# Printing out info (no BED file)
def printInfo_minimal(options):
    print 'Input, output and settings:'
    print "--------------------------------------------------------------------------------------"
    print "Input file name:        " + options.input
    print ''
    print "Output formats:         _summary"
    print "Output files prefix:    " + options.output
    print "--------------------------------------------------------------------------------------"

# Setting default values for unspecified configuration options
def defaultConfigs(config):
    if not 'duplicates' in config.keys(): config['duplicates'] = True
    if not 'outputs' in config.keys(): config['outputs'] = {'regions': True, 'profiles': True, 'gui': False}
    if not 'regions' in config['outputs'].keys(): config['outputs']['regions'] = True
    if not 'profiles' in config['outputs'].keys(): config['outputs']['profiles'] = True
    if not 'gui' in config['outputs'].keys(): config['outputs']['gui'] = False
    if not 'low_bq' in config.keys(): config['low_bq'] = 10
    if not 'low_mq' in config.keys(): config['low_mq'] = 20
    if not 'only_fail_profiles' in config.keys(): config['only_fail_profiles'] = False
    if not 'transcript_db' in config.keys(): config['transcript_db'] = None
    if not 'pass' in config.keys() or len(config['pass']) == 0: config['pass'] = None
    if not 'transcript' in config.keys(): config['transcript'] = {'regions': True, 'profiles': True, 'poor': True}
    if not 'regions' in config['transcript'].keys(): config['transcript']['regions'] = True
    if not 'profiles' in config['transcript'].keys(): config['transcript']['profiles'] = True
    if not 'poor' in config['transcript'].keys(): config['transcript']['poor'] = True
    if not 'direction' in config.keys(): config['direction'] = False

# Creating target name list
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

# Finding break points in the input file
def findFileBreaks(inputf, threads):
    ret = []
    started = False
    counter = 0
    for line in open(inputf):
        counter += 1
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        if not started:
            started = True
            first = counter
    delta = int((counter - first + 1) / threads)

    counter = 0
    prevchrom = ''
    i = 0
    start = 0
    for line in open(inputf):
        counter += 1
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue

        cols = line.split('\t')
        chrom = cols[0]

        if counter > delta and not chrom == prevchrom:
            ret.append((start, i))
            start = i + 1
            counter = 0

        prevchrom = chrom
        i += 1
    ret.append((start, ''))

    return ret

# Merging tmp files to final output files
def mergeTmpFiles(options, config):

    # Merging _region output tmp files
    if config['outputs']['regions']:

        filenames = []
        for i in range(1, int(options.threads) + 1):
            filenames.append(options.output + '_regions_tmp_' + str(i) + '.txt')

        with open(options.output + '_regions.txt', 'w') as outfile:
            targetheader = ['Region', 'Chromosome', 'Start_position', 'End_position']
            if config['transcript']['regions'] and not config['transcript_db'] is None: targetheader.extend(
                ['Start_transcript', 'End_transcript'])
            if not config['pass'] is None: targetheader.append('Pass_or_fail')
            targetheader.extend(['RC', 'MEDCOV', 'MINCOV', 'MEDQCOV', 'MINQCOV', 'MAXFLMQ', 'MAXFLBQ'])

            if config['direction']:
                targetheader.extend(['MEDCOV+', 'MINCOV+', 'MEDQCOV+', 'MINQCOV+', 'MAXFLMQ+', 'MAXFLBQ+'])
                targetheader.extend(['MEDCOV-', 'MINCOV-', 'MEDQCOV-', 'MINQCOV-', 'MAXFLMQ-', 'MAXFLBQ-'])

            outfile.write('#' + '\t'.join(targetheader) + '\n')

            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        for fn in filenames:
            os.remove(fn)

    # Merging _profiles output tmp files
    if config['outputs']['profiles']:

        filenames = []
        for i in range(1, int(options.threads) + 1):
            filenames.append(options.output + '_profiles_tmp_' + str(i) + '.txt')

        with open(options.output + '_profiles.txt', 'w') as outfile:
            profheader = ['Chromosome', 'Position']
            if config['transcript']['profiles'] and not config['transcript_db'] is None: profheader.append(
                'Transcript_coordinate')
            profheader.extend(['COV', 'QCOV', 'MEDBQ', 'FLBQ', 'MEDMQ', 'FLMQ'])

            if config['direction']:
                profheader.extend(['COV+', 'QCOV+', 'MEDBQ+', 'FLBQ+', 'MEDMQ+', 'FLMQ+'])
                profheader.extend(['COV-', 'QCOV-', 'MEDBQ-', 'FLBQ-', 'MEDMQ-', 'FLMQ-'])

            outfile.write('#' + '\t'.join(profheader) + '\n')

            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        for fn in filenames: os.remove(fn)

    # Merging gui data json tmp files
    if config['outputs']['gui']:

        filenames = []
        for i in range(1, int(options.threads) + 1):
            filenames.append(options.output + '_results_tmp_' + str(i) + '.json')

        with open(options.output + '_gui/data/results.js', 'w') as outfile:
            outfile.write('function readData() {\n')
            outfile.write('\tdata={\"targets\":[')

            for i in range(len(filenames)):
                fname = filenames[i]
                if i > 0: outfile.write(',')
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

            outfile.write(']')

        for fn in filenames: os.remove(fn)

    # Merging _poor output tmp files
    if config['transcript']['poor'] and not config['transcript_db'] is None and config['outputs']['profiles']:

        filenames = []
        for i in range(1, int(options.threads) + 1):
            filenames.append(options.output + '_poor_tmp_' + str(i) + '.txt')

        with open(options.output + '_poor.txt', 'w') as outfile:
            poorheader = ['Region', 'Chromosome', 'Start_position', 'End_position', 'Start_transcript',
                          'End_transcript']
            outfile.write('#' + '\t'.join(poorheader) + '\n')

            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        for fn in filenames: os.remove(fn)

# Writing _summary output file
def output_summary(options, chromdata):
    out_summary = open(options.output + '_summary.txt', 'w')
    out_summary.write('#CHROM\tRC\tRCIN\tRCOUT\n')
    mapped = chromdata['Mapped']
    unmapped = chromdata['Unmapped']
    total = chromdata['Total']
    out_summary.write('Total\t' + str(total) + '\t-\t-\n')
    out_summary.write('Unmapped\t' + str(unmapped) + '\t-\t-\n')
    record = 'Mapped' + '\t' + str(mapped['RC']) + '\t' + str(mapped['RCIN']) + '\t' + str(mapped['RCOUT'])
    out_summary.write(record + '\n')
    for i in range(len(chromdata['Chroms'])):
        chromres = chromdata['Chroms'][i]
        chrom = chromres['CHROM']
        record = chrom + '\t' + str(chromres['RC']) + '\t' + str(chromres['RCIN']) + '\t' + str(chromres['RCOUT'])
        out_summary.write(record + '\n')
    out_summary.close()

# Calculating chromosome data
def calculateChromdata(samfile, ontarget):
    # Initializing progress info
    sys.stdout.write('\rFinalizing analysis ... 0.0%')
    sys.stdout.flush()

    chromdata = dict()
    chroms = samfile.references
    chromsres = []
    alltotal = 0
    allon = 0
    alloff = 0

    # Iterate through chromosomes
    i = 0
    for chrom in chroms:

        total = sum(1 for _ in samfile.fetch(chrom))

        if 'chr' + chrom in ontarget.keys():
            on = int(ontarget['chr' + chrom])
            off = total - on
        else:
            on = 0
            off = 0
        chromsres.append({'CHROM': chrom, 'RC': total, 'RCIN': on, 'RCOUT': off})
        alltotal += total
        allon += on
        alloff += off

        i += 1

        # Updating progress info
        x = round(100 * i / len(chroms), 1)
        x = min(x, 100.0)
        sys.stdout.write('\rFinalizing analysis ... ' + str(x) + '%')
        sys.stdout.flush()

    chromdata['Chroms'] = chromsres
    chromdata['Mapped'] = {'RC': alltotal, 'RCIN': allon, 'RCOUT': alloff}

    allreads = pysam.flagstat(options.input)[0]
    allreads = allreads[:allreads.find('+')]
    allreads = int(allreads.strip())

    chromdata['Total'] = allreads
    chromdata['Unmapped'] = allreads - alltotal

    # Finalizing progress info
    sys.stdout.write('\rFinalizing analysis ... 100.0% - Done')
    sys.stdout.flush()
    print ''

    return chromdata

# Calculating chromosome data (minimal mode)
def calculateChromdata_minimal(samfile):
    print ''
    # Initializing progress info
    sys.stdout.write('\rRunning analysis ... 0.0%')
    sys.stdout.flush()

    chromdata = dict()
    chroms = samfile.references
    chromsres = []
    alltotal = 0

    i = 0
    # Iterate through chromosomes
    for chrom in chroms:
        total = sum(1 for _ in samfile.fetch(chrom))
        chromsres.append({'CHROM': chrom, 'RC': total})
        alltotal += total

        i += 1

        # Updating progress info
        x = round(100 * i / len(chroms), 1)
        x = min(x, 100.0)
        sys.stdout.write('\rRunning analysis ... ' + str(x) + '%')
        sys.stdout.flush()

    chromdata['Chroms'] = chromsres
    chromdata['Mapped'] = {'RC': alltotal}

    allreads = pysam.flagstat(options.input)[0]
    allreads = allreads[:allreads.find('+')]
    allreads = int(allreads.strip())

    chromdata['Total'] = allreads
    chromdata['Unmapped'] = allreads - alltotal

    # Finalizing progress info
    sys.stdout.write('\rRunning analysis ... 100.0% - Done')
    sys.stdout.flush()
    print ''

    return chromdata

# Writing _summary output file (minimal mode)
def output_summary_minimal(options, chromdata):
    out_summary = open(options.output + '_summary.txt', 'w')
    out_summary.write('#CHROM\tRC\n')
    mapped = chromdata['Mapped']
    unmapped = chromdata['Unmapped']
    total = chromdata['Total']
    out_summary.write('Total\t' + str(total) + '\n')
    out_summary.write('Unmapped\t' + str(unmapped) + '\n')
    record = 'Mapped' + '\t' + str(mapped['RC'])
    out_summary.write(record + '\n')
    for i in range(len(chromdata['Chroms'])):
        chromres = chromdata['Chroms'][i]
        chrom = chromres['CHROM']
        record = chrom + '\t' + str(chromres['RC'])
        out_summary.write(record + '\n')
    out_summary.close()


def finalizeJSONOutput(options):

    out_json = open(options.output + '_gui/data/results.js', 'a')

    newchromsres = []
    others = {'CHROM': '...', 'RC': 0, 'RCIN': 0, 'RCOUT': 0}
    chrnames = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                '20', '21', '22', 'X', 'Y', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                'chr21', 'chr22', 'chrX', 'chrY']
    chromsres = chromdata['Chroms']
    for x in chromsres:
        if x['CHROM'] in chrnames:
            newchromsres.append(x)
        else:
            others['RC'] += x['RC']
            others['RCIN'] += x['RCIN']
            others['RCOUT'] += x['RCOUT']
    newchromsres.append(others)
    chromdata['Chroms'] = newchromsres

    out_json.write(',\"chromdata\":' + json.dumps(chromdata, separators=(',', ':')))

    infn = options.input
    if '/' in infn:
        infn = infn[infn.rfind('/') + 1:]
    if len(infn) > 27: infn = infn[:27] + '...bam'
    out_json.write(',\"input\":\"' + infn + '\"')

    out_json.write(',\"direction\":' + str(config['direction']).lower())
    out_json.write(',\"duplicates\":' + str(config['duplicates']).lower())
    out_json.write(',\"ntargets\":' + str(numOfTargets))
    out_json.write(',\"unique\":' + str(uniqueIDs))
    out_json.write(',\"nfailed\":' + str(failedtargets))
    out_json.write(',\"uniquefailed\":' + str(len(uniqueids)))

    bedfn = options.bedfile
    if '/' in bedfn:
        bedfn = bedfn[bedfn.rfind('/') + 1:]
    if len(bedfn) > 20: bedfn = bedfn[:20] + '...bed'
    out_json.write(',\"bedfile\":\"' + bedfn + '\"')

    now = datetime.datetime.now()
    out_json.write(',\"date\":\"' + now.strftime("%d-%m-%Y, %H:%M") + "\"")

    passdef = []
    for k, v in config['pass'].iteritems():
        [met,minmax] = k.split('_')
        if minmax == 'MIN':
            passdef.append(met + '>' + str(v))
        else:
            passdef.append(met + '<' + str(v))
    passdefstr = ', '.join(passdef)
    out_json.write(',\"passdef\":\"' + passdefstr + '\"')

    passmets = dict()
    for k, v in config['pass'].iteritems():
        [met, minmax] = k.split('_')
        if met.endswith('QCOV'):  passmets['QCOV'] = v
        elif met.endswith('COV'):  passmets['COV'] = v
        elif met.endswith('FLMQ'):  passmets['FLMQ'] = v
        elif met.endswith('FLBQ'): passmets['FLBQ'] = v
        passmets[met] = v
    out_json.write(',\"passmets\":' + json.dumps(passmets, separators=(',', ':')))

    out_json.write('}\n')
    out_json.write('\treturn data\n')
    out_json.write('}\n')
    out_json.close()

#########################################################################################################################################


numpy.seterr(all='ignore')
warnings.simplefilter("ignore", RuntimeWarning)

print ""
print "======================================================================================"
print 'CoverView v1.2.0 started running: ', datetime.datetime.now()
print ""

# Command line argument parsing
parser = OptionParser(usage='usage: python %prog [options]')
parser.add_option("-i", "--input", default='input.bam', dest='input', action='store', help="Input (BAM) filename [default value: %default]")
parser.add_option("-o", "--output", default='output', dest='output', action='store', help="Output filename [default value: %default]")
parser.add_option("-b", "--bed", default=None, dest='bedfile', action='store', help="Input BED filename [default value: %default]")
parser.add_option("-c", "--config", default=None, dest='config', action='store', help="Configuration file [default value: %default]")
parser.add_option("-t", "--threads", default=1, dest='threads', action='store', help="Number of processes used [default value: %default]")
(options, args) = parser.parse_args()

# Loading configuration file
config = dict()
if not options.config is None:
    with open(options.config) as config_file: config = json.load(config_file)
defaultConfigs(config)

# Removing unremoved temporary files if exist
if os.path.isfile(options.output + '_failedtargets.txt'): os.remove(options.output + '_failedtargets.txt')
if os.path.isfile(options.output + '_reads_on_target.txt'): os.remove(options.output + '_reads_on_target.txt')

# Minimal mode (no BED file)
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

# Creating output directory for gui
if config['outputs']['gui']:
    dir = options.output + '_gui'
    if os.path.exists(dir): shutil.rmtree(dir)
    os.makedirs(dir)
    os.makedirs(dir+'/data')
    cvdir = os.path.dirname(os.path.realpath(__file__))
    shutil.copy(cvdir+'/gui.html', dir+'/'+options.output+'_coverview.html')
    shutil.copytree(cvdir+'/lib', dir+'/lib')

# Creating target name list
names, uniqueIDs = makeNames(options.bedfile)
numOfTargets = len(names)

# Setting sample name
samplename = options.input[:options.input.rfind('.bam')]

# Defining break points in the input BED file
breaks = findFileBreaks(options.bedfile, int(options.threads))
options.threads = len(breaks)

# Print out info
printInfo(options, config, numOfTargets)

# Initializing processes
threadidx = 0
processes = []
for (startline, endline) in breaks:
    threadidx += 1
    processes.append(SingleJob(threadidx, options, config, startline, endline, names))

# Running processes
for process in processes: process.start()
for process in processes: process.join()

# Merging tmp files
if int(options.threads) > 1: mergeTmpFiles(options, config)

# Reading on-target read counts from tmp file
ontarget = dict()
for threadidx in range(1, int(options.threads) + 1):
    for line in open(options.output + '_reads_on_target_' + str(threadidx) + '.txt'):
        [key, value] = line.split(':')
        if key in ontarget.keys(): ontarget[key] += int(value.strip())
        else: ontarget[key] = int(value.strip())
    os.remove(options.output + '_reads_on_target_' + str(threadidx) + '.txt')

# Reading number of failed targets from tmp file
failedtargets = 0
uniqueids = set()
for threadidx in range(1, int(options.threads) + 1):
    for line in open(options.output + '_failedtargets_' + str(threadidx) + '.txt'):
        line = line.strip()
        fails = line[1:-1].split(',')
        failedtargets += len(fails)
        for x in fails:
            x = x.strip()
            uniqueids.add(x[1:-1])
    os.remove(options.output + '_failedtargets_' + str(threadidx) + '.txt')

# Closing progress info
if failedtargets > 0: print ' - Done. (' + str(failedtargets) + ' failed regions)'
else: print ' - Done.'

# Calculate and output chromosome summary data
samfile = pysam.Samfile(options.input, "rb")
chromdata = calculateChromdata(samfile, ontarget)
output_summary(options, chromdata)

# Finalizing JSON output file
if config['outputs']['gui']: finalizeJSONOutput(options)

# Goodbye message
print ""
print 'CoverView v1.2.0 succesfully finished: ', datetime.datetime.now()
print "======================================================================================"
print ""

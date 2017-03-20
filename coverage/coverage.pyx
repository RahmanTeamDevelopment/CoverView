"""
Utilities for calculating coverage summaries from a BAM file
"""

from __future__ import division

import numpy
from statistics cimport QualityHistogram

class SimpleCoverageCalculator(object):
    """
    Utility class for computing coverage summaries for the simple case when we do
    not care about direction.
    """
    def __init__(self, chrom, begin, end, bq_cutoff, mq_cutoff):
        self.begin = begin
        self.end = end
        self.bq_cutoff = bq_cutoff
        self.mq_cutoff = mq_cutoff
        self.ret_COV = [0] * (end - begin)
        self.ret_QCOV = [0] * (end - begin)
        self.ret_MEDBQ = [float('NaN')] * (end - begin)
        self.ret_FLBQ = [float('NaN')] * (end - begin)
        self.ret_MEDMQ = [float('NaN')] * (end - begin)
        self.ret_FLMQ = [float('NaN')] * (end - begin)

        self.bq_hists = []
        self.mq_hists = []

        for x in self.ret_COV:
            self.bq_hists.append(QualityHistogram())
            self.mq_hists.append(QualityHistogram())

    def add_reads(self, reads):
        """
        """
        # Caching object variables to local variables here, as these are used
        # very frequently in this function.
        cdef QualityHistogram bq_hist
        cdef QualityHistogram mq_hist

        region_size = self.end - self.begin
        begin = self.begin
        end = self.end
        bq_hists = self.bq_hists
        mq_hists = self.mq_hists
        bq_cutoff = self.bq_cutoff
        mq_cutoff = self.mq_cutoff
        ret_COV = self.ret_COV
        ret_QCOV = self.ret_QCOV

        for read in reads:
            mapping_quality = read.mapping_quality
            base_qualities = read.query_qualities

            for index, ref_pos in enumerate(read.get_reference_positions(full_length=True)):

                if ref_pos is None:
                    continue

                # TODO: add 1 to cov for deletions and 1 to qco for high map qual deletions
                if begin <= ref_pos < end:
                    offset = ref_pos - begin
                    base_quality = base_qualities[index]
                    bq_hists[offset].add_data(base_quality)
                    mq_hists[offset].add_data(mapping_quality)
                    ret_COV[offset] += 1

                    if mapping_quality >= mq_cutoff and base_quality >= bq_cutoff:
                        ret_QCOV[offset] += 1

        for i in xrange(len(self.bq_hists)):

            bq_hist = self.bq_hists[i]
            mq_hist = self.mq_hists[i]

            self.ret_MEDBQ[i] = bq_hist.compute_median()
            self.ret_MEDMQ[i] = mq_hist.compute_median()
            self.ret_FLBQ[i] = round(bq_hist.compute_fraction_below_threshold(int(self.bq_cutoff)), 3)
            self.ret_FLMQ[i] = round(mq_hist.compute_fraction_below_threshold(int(self.mq_cutoff)), 3)

    def add_pileup(self, pileup_column, i):
        self.ret_COV[i] = pileup_column.n
        bqs = []
        mqs = []
        qcov = 0

        for pileupread in pileup_column.pileups:
            if pileupread.is_del:
                if pileupread.alignment.mapq >= self.mq_cutoff:
                    qcov += 1
                continue

            bq = pileupread.alignment.query_qualities[pileupread.query_position]
            bqs.append(bq)
            mqs.append(pileupread.alignment.mapq)

            if pileupread.alignment.mapq >= self.mq_cutoff and bq >= self.bq_cutoff:
                qcov += 1

        self.ret_QCOV[i] = qcov
        self.ret_MEDBQ[i] = numpy.median(bqs)
        self.ret_MEDMQ[i] = numpy.median(mqs)

        if len(bqs) > 0:
            self.ret_FLBQ[i] = round(len([x for x in bqs if x < self.bq_cutoff]) / len(bqs), 3)

        if len(mqs) > 0:
            self.ret_FLMQ[i] = round(len([x for x in mqs if x < self.mq_cutoff]) / len(mqs), 3)

    def get_coverage_summary(self):
        return (
            self.ret_COV,
            self.ret_QCOV,
            self.ret_MEDBQ,
            self.ret_FLBQ,
            self.ret_MEDMQ,
            self.ret_FLMQ
        )


class DirectionalCoverageCalculator(object):
    """
    Utility class for computing coverage summaries for the case when we do
    care about direction.
    """
    def __init__(self, chrom, begin, end, bq_cutoff, mq_cutoff):
        self.bq_cutoff = bq_cutoff
        self.mq_cutoff = mq_cutoff
        self.ret_COV = [0] * (end - begin)
        self.ret_QCOV = [0] * (end - begin)
        self.ret_MEDBQ = [float('NaN')] * (end - begin)
        self.ret_FLBQ = [float('NaN')] * (end - begin)
        self.ret_MEDMQ = [float('NaN')] * (end - begin)
        self.ret_FLMQ = [float('NaN')] * (end - begin)
        self.ret_COV_f = [0] * (end - begin)
        self.ret_QCOV_f = [0] * (end - begin)
        self.ret_MEDBQ_f = [float('NaN')] * (end - begin)
        self.ret_FLBQ_f = [float('NaN')] * (end - begin)
        self.ret_MEDMQ_f = [float('NaN')] * (end - begin)
        self.ret_FLMQ_f = [float('NaN')] * (end - begin)
        self.ret_COV_r = [0] * (end - begin)
        self.ret_QCOV_r = [0] * (end - begin)
        self.ret_MEDBQ_r = [float('NaN')] * (end - begin)
        self.ret_FLBQ_r = [float('NaN')] * (end - begin)
        self.ret_MEDMQ_r = [float('NaN')] * (end - begin)
        self.ret_FLMQ_r = [float('NaN')] * (end - begin)

    def add_pileup(self, pileup_column, i):
        self.ret_COV[i] = pileup_column.n
        bqs = []
        mqs = []
        qcov = 0
        cov_f = 0
        cov_r = 0
        bqs_f = []
        bqs_r = []
        mqs_f = []
        mqs_r = []
        qcov_f = 0
        qcov_r = 0

        for pileupread in pileup_column.pileups:

            if pileupread.alignment.is_reverse:
                cov_r += 1
            else:
                cov_f += 1

            if pileupread.is_del:
                if pileupread.alignment.mapq >= self.mq_cutoff:
                    qcov += 1
                    if pileupread.alignment.is_reverse:
                        qcov_r += 1
                    else:
                        qcov_f += 1
                continue

            bq = pileupread.alignment.query_qualities[pileupread.query_position]
            bqs.append(bq)

            if pileupread.alignment.is_reverse:
                bqs_r.append(bq)
            else:
                bqs_f.append(bq)

            mqs.append(pileupread.alignment.mapq)

            if pileupread.alignment.is_reverse:
                mqs_r.append(pileupread.alignment.mapq)
            else:
                mqs_f.append(pileupread.alignment.mapq)

            if pileupread.alignment.mapq >= self.mq_cutoff and bq >= self.bq_cutoff:
                qcov += 1
                if pileupread.alignment.is_reverse:
                    qcov_r += 1
                else:
                    qcov_f += 1

        self.ret_QCOV[i] = qcov
        self.ret_MEDBQ[i] = numpy.median(bqs)
        self.ret_MEDMQ[i] = numpy.median(mqs)

        if len(bqs) > 0:
            self.ret_FLBQ[i] = round(len([x for x in bqs if x < self.bq_cutoff]) / len(bqs), 3)

        if len(mqs) > 0:
            self.ret_FLMQ[i] = round(len([x for x in mqs if x < self.mq_cutoff]) / len(mqs), 3)

        self.ret_COV_f[i] = cov_f
        self.ret_COV_r[i] = cov_r
        self.ret_QCOV_f[i] = qcov_f
        self.ret_QCOV_r[i] = qcov_r
        self.ret_MEDBQ_f[i] = numpy.median(bqs_f)
        self.ret_MEDBQ_r[i] = numpy.median(bqs_r)
        self.ret_MEDMQ_f[i] = numpy.median(mqs_f)
        self.ret_MEDMQ_r[i] = numpy.median(mqs_r)

        if len(bqs_f) > 0:
            self.ret_FLBQ_f[i] = round(len([x for x in bqs_f if x < self.bq_cutoff]) / len(bqs_f), 3)

        if len(bqs_r) > 0:
            self.ret_FLBQ_r[i] = round(len([x for x in bqs_r if x < self.bq_cutoff]) / len(bqs_r), 3)

        if len(mqs_f) > 0:
            self.ret_FLMQ_f[i] = round(len([x for x in mqs_f if x < self.mq_cutoff]) / len(mqs_f), 3)

        if len(mqs_r) > 0:
            self.ret_FLMQ_r[i] = round(len([x for x in mqs_r if x < self.mq_cutoff]) / len(mqs_r), 3)

    def get_coverage_summary(self):
        return (
            self.ret_COV,
            self.ret_QCOV,
            self.ret_MEDBQ,
            self.ret_FLBQ,
            self.ret_MEDMQ,
            self.ret_FLMQ,
            self.ret_COV_f,
            self.ret_QCOV_f,
            self.ret_MEDBQ_f,
            self.ret_FLBQ_f,
            self.ret_MEDMQ_f,
            self.ret_FLMQ_f,
            self.ret_COV_r,
            self.ret_QCOV_r,
            self.ret_MEDBQ_r,
            self.ret_FLBQ_r,
            self.ret_MEDMQ_r,
            self.ret_FLMQ_r
        )


def get_valid_chromosome_name(chrom, bam_file):
    """
    Return a chromosome name that matches the naming convention used in the BAM file, i.e.
    add or remove 'chr' as appropriate.
    """
    chr_prefix = bam_file.references[0].startswith('chr')

    if chr_prefix and not chrom.startswith('chr'):
        return 'chr' + chrom

    if not chr_prefix and chrom.startswith('chr'):
        return chrom[3:]

    return chrom


def get_profiles(bam_file, region, config):
    """
    Calculate and return coverage metrics for a specified region. Metrics include total
    coverage, coverage above the required base-quality and mapping quality threshold, fractions of
    low base qualities at a given position, fractions of low mapping quality reads covering a given
    position, median mapping qualities and median base qualities.

    This is by far the most computationally expensive part of CoverView. > 90% of the run-time is
    currently spent in this function.
    """
    bq_cutoff = float(config['low_bq'])
    mq_cutoff = float(config['low_mq'])

    chrom, beg_end = region.split(":")
    begin = int(beg_end.split("-")[0]) - 1
    end = int(beg_end.split("-")[1]) - 1

    good_chrom = get_valid_chromosome_name(chrom, bam_file)

    if not good_chrom in bam_file.references:
        return None

    coverage_calc = SimpleCoverageCalculator(good_chrom, begin, end, bq_cutoff, mq_cutoff)
    reads = bam_file.fetch(good_chrom, begin, end)
    coverage_calc.add_reads(reads)
    # if config['direction']:
    #     coverage_calc = DirectionalCoverageCalculator(chrom, begin, end, bq_cutoff, mq_cutoff)
    # else:
    #     coverage_calc = SimpleCoverageCalculator(chrom, begin, end, bq_cutoff, mq_cutoff)
    #
    # if config['duplicates']:
    #     x = bam_file.pileup(good_chrom, begin + 1, end + 1, mask=0, truncate=True)
    # else:
    #     x = bam_file.pileup(good_chrom, begin + 1, end + 1, truncate=True)
    #
    # for i, pileup_column in enumerate(x):
    #     coverage_calc.add_pileup(pileup_column, i)

    return coverage_calc.get_coverage_summary()
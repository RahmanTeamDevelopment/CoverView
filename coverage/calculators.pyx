"""
Utilities for calculating coverage summaries from a BAM file
"""

from __future__ import division

import numpy
from coverage.statistics cimport QualityHistogram

from libc.stdint cimport uint32_t, uint64_t
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment, bam1_t, BAM_CIGAR_MASK,\
    BAM_CIGAR_SHIFT, BAM_CINS, BAM_CSOFT_CLIP, BAM_CREF_SKIP, BAM_CMATCH, BAM_CDEL

from cpython cimport array

cdef extern from "hts.h":
    ctypedef struct hts_idx_t
    int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped)


cdef class SimpleCoverageCalculator(object):
    """
    Utility class for computing coverage summaries for the simple case when we do
    not care about direction.
    """
    cdef int begin
    cdef int end
    cdef int bq_cutoff
    cdef int mq_cutoff
    cdef array.array COV,QCOV,MEDBQ,FLBQ,MEDMQ,FLMQ
    cdef list bq_hists
    cdef list mq_hists

    def __init__(self, chrom, begin, end, bq_cutoff, mq_cutoff):
        self.begin = begin
        self.end = end
        self.bq_cutoff = bq_cutoff
        self.mq_cutoff = mq_cutoff
        self.COV = array.array('l', [0] * (end - begin))
        self.QCOV = array.array('l', [0] * (end - begin))
        self.MEDBQ = array.array('f', [float('NaN')] * (end - begin))
        self.FLBQ = array.array('f', [float('NaN')] * (end - begin))
        self.MEDMQ = array.array('f', [float('NaN')] * (end - begin))
        self.FLMQ = array.array('f', [float('NaN')] * (end - begin))

        self.bq_hists = []
        self.mq_hists = []

        for x in self.COV:
            self.bq_hists.append(QualityHistogram())
            self.mq_hists.append(QualityHistogram())

    def add_reads(self, reads):
        """
        """
        # Caching object variables to local variables here, as these are used
        # very frequently in this function.
        cdef QualityHistogram bq_hist
        cdef QualityHistogram mq_hist
        cdef int region_size = self.end - self.begin
        cdef int begin = self.begin
        cdef int end = self.end
        cdef list bq_hists = self.bq_hists
        cdef list mq_hists = self.mq_hists
        cdef int bq_cutoff = self.bq_cutoff
        cdef int mq_cutoff = self.mq_cutoff
        cdef long* COV = self.COV.data.as_longs
        cdef long* QCOV = self.QCOV.data.as_longs
        cdef AlignedSegment read
        cdef int index
        cdef int base_quality
        cdef int mapping_quality
        cdef array.array base_qualities_array
        cdef unsigned char* base_qualities
        cdef uint32_t k, i, pos, n_cigar
        cdef int op,l
        cdef uint32_t* cigar_p
        cdef bam1_t* src

        for read in reads:
            mapping_quality = read.mapping_quality
            base_qualities_array = read.query_qualities
            base_qualities = base_qualities_array.data.as_uchars
            src = read._delegate
            n_cigar = src.core.n_cigar

            if n_cigar == 0:
                break

            pos = src.core.pos
            cigar_p = <uint32_t*> (src.data + src.core.l_qname)
            index = 0

            for k from 0 <= k < n_cigar:
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT

                if op == BAM_CSOFT_CLIP or op == BAM_CINS:
                    index += l
                elif op == BAM_CMATCH:
                    for i from pos <= i < pos + l:
                        # TODO: add 1 to cov for deletions and 1 to qco for high map qual deletions
                        if begin < i <= end:
                            offset = i - (begin + 1)
                            base_quality = base_qualities[index + (i-pos)]
                            bq_hist = bq_hists[offset]
                            mq_hist = mq_hists[offset]
                            bq_hist.add_data(base_quality)
                            mq_hist.add_data(mapping_quality)
                            COV[offset] += 1

                            if mapping_quality >= mq_cutoff and base_quality >= bq_cutoff:
                                QCOV[offset] += 1

                    pos += l
                    index += l
                elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                    pos += l

        for i in xrange(len(self.bq_hists)):
            bq_hist = self.bq_hists[i]
            mq_hist = self.mq_hists[i]
            self.MEDBQ[i] = bq_hist.compute_median()
            self.MEDMQ[i] = mq_hist.compute_median()
            self.FLBQ[i] = round(bq_hist.compute_fraction_below_threshold(int(self.bq_cutoff)), 3)
            self.FLMQ[i] = round(mq_hist.compute_fraction_below_threshold(int(self.mq_cutoff)), 3)

    def add_pileup(self, pileup_column, i):
        self.COV[i] = pileup_column.n
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

        self.QCOV[i] = qcov
        self.MEDBQ[i] = numpy.median(bqs)
        self.MEDMQ[i] = numpy.median(mqs)

        if len(bqs) > 0:
            self.FLBQ[i] = round(len([x for x in bqs if x < self.bq_cutoff]) / len(bqs), 3)

        if len(mqs) > 0:
            self.FLMQ[i] = round(len([x for x in mqs if x < self.mq_cutoff]) / len(mqs), 3)

    def get_coverage_summary(self):
        return (
            list(self.COV),
            list(self.QCOV),
            list(self.MEDBQ),
            list(self.FLBQ),
            list(self.MEDMQ),
            list(self.FLMQ)
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


def get_num_mapped_reads_covering_chromosome(bam_file, chrom):
    """
    Return the total number of reads covering the specified chromsome
    in the specified BAM file. This is an optimisation, which makes use of the
    index statistics in the BAM index, which record the nunmber of alignments. This is
    an O(1) operation rather than the O(N) operation of looping through all reads in the
    file with read.tid == chromID.
    """
    cdef AlignmentFile the_file = bam_file
    cdef hts_idx_t* index = the_file.index
    cdef uint64_t n_mapped = 0
    cdef uint64_t n_unmapped = 0
    cdef int tid = the_file.get_tid(chrom)

    hts_idx_get_stat(index, tid, &n_mapped, &n_unmapped)
    return n_mapped


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
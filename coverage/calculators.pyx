"""
Utilities for calculating coverage summaries from a BAM file
"""

from __future__ import division

import numpy
from coverage.statistics cimport QualityHistogramArray

from libc.stdint cimport uint32_t, uint64_t, uint8_t

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment, bam1_t, BAM_CIGAR_MASK,\
    BAM_CIGAR_SHIFT, BAM_CINS, BAM_CSOFT_CLIP, BAM_CREF_SKIP, BAM_CMATCH, BAM_CDEL, BAM_FDUP

from pysam.libcalignmentfile cimport IteratorRowRegion
from cpython cimport array


cdef extern from "hts.h":
    ctypedef struct hts_idx_t
    ctypedef struct hts_itr_t
    ctypedef struct BGZF
    ctypedef struct htsFile
    int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped)
    BGZF *hts_get_bgzfp(htsFile *fp)
    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)


cdef extern from "sam.h":
    uint8_t* bam_get_qual(bam1_t* b)


cdef class SimpleCoverageCalculator(object):
    """
    Utility class for computing coverage summaries for the simple case when we do
    not care about direction.
    """
    cdef int begin
    cdef int end
    cdef int bq_cutoff
    cdef int mq_cutoff
    cdef int n_reads_in_region
    cdef array.array COV,QCOV,MEDBQ,FLBQ,MEDMQ,FLMQ
    cdef QualityHistogramArray bq_hists
    cdef QualityHistogramArray mq_hists

    def __init__(self, chrom, begin, end, bq_cutoff, mq_cutoff):
        cdef int bases_in_region = end - begin
        self.begin = begin
        self.end = end
        self.bq_cutoff = bq_cutoff
        self.mq_cutoff = mq_cutoff
        self.n_reads_in_region = 0
        self.COV = array.array('l', [0] * (bases_in_region))
        self.QCOV = array.array('l', [0] * (bases_in_region))
        self.MEDBQ = array.array('f', [float('NaN')] * (bases_in_region))
        self.FLBQ = array.array('f', [float('NaN')] * (bases_in_region))
        self.MEDMQ = array.array('f', [float('NaN')] * (bases_in_region))
        self.FLMQ = array.array('f', [float('NaN')] * (bases_in_region))
        self.bq_hists = QualityHistogramArray(bases_in_region)
        self.mq_hists = QualityHistogramArray(bases_in_region)

    cdef void add_reads(self, IteratorRowRegion read_iterator):
        """
        """
        cdef int region_size = self.end - self.begin
        cdef int begin = self.begin
        cdef int end = self.end
        cdef int bq_cutoff = self.bq_cutoff
        cdef int mq_cutoff = self.mq_cutoff
        cdef long* COV = self.COV.data.as_longs
        cdef long* QCOV = self.QCOV.data.as_longs
        cdef int index
        cdef int base_quality
        cdef int mapping_quality
        cdef uint32_t k, i
        cdef uint32_t pos, n_cigar
        cdef int op,l
        cdef uint32_t* cigar_p
        cdef bam1_t* src
        cdef uint8_t* base_qualities
        cdef int iterator_status = 0

        while True:
            iterator_status = hts_itr_next(
                hts_get_bgzfp(read_iterator.htsfile),
                read_iterator.iter,
                read_iterator.b,
                read_iterator.htsfile
            )

            if iterator_status < 0:
                break

            src = read_iterator.b

            # Skip duplicate reads
            #if src.core.flag & BAM_FDUP != 0:
            #    continue

            n_cigar = src.core.n_cigar

            if n_cigar == 0:
                break

            self.n_reads_in_region += 1
            mapping_quality = src.core.qual
            base_qualities = bam_get_qual(src)
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
                            self.bq_hists.add_data(offset, base_quality)
                            self.mq_hists.add_data(offset, mapping_quality)
                            COV[offset] += 1

                            if mapping_quality >= mq_cutoff and base_quality >= bq_cutoff:
                                QCOV[offset] += 1

                    pos += l
                    index += l
                elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                    for i from pos <= i < pos + l:
                        if begin < i <= end:
                            offset = i - (begin + 1)
                            COV[offset] += 1
                            if mapping_quality >= mq_cutoff:
                                QCOV[offset] += 1
                    pos += l

        self.compute_summary_statistics_for_region()

    cdef void compute_summary_statistics_for_region(self):
        cdef float* MEDBQ = self.MEDBQ.data.as_floats
        cdef float* MEDMQ = self.MEDMQ.data.as_floats
        cdef float* FLBQ = self.FLBQ.data.as_floats
        cdef float* FLMQ = self.FLMQ.data.as_floats

        for i in xrange(self.end - self.begin):
            MEDBQ[i] = self.bq_hists.compute_median(i)
            MEDMQ[i] = self.mq_hists.compute_median(i)
            FLBQ[i] = self.bq_hists.compute_fraction_below_threshold(i, <int>(self.bq_cutoff))
            FLMQ[i] = self.mq_hists.compute_fraction_below_threshold(i, <int>(self.mq_cutoff))

    def get_coverage_summary(self):
        return (
            self.n_reads_in_region,
            self.COV,
            self.QCOV,
            self.MEDBQ,
            self.FLBQ,
            self.MEDMQ,
            self.FLMQ
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
    end = int(beg_end.split("-")[1]) - 1 + 1 # TODO: sort out the indexing!

    good_chrom = get_valid_chromosome_name(chrom, bam_file)

    if not good_chrom in bam_file.references:
        return None

    coverage_calc = SimpleCoverageCalculator(good_chrom, begin, end, bq_cutoff, mq_cutoff)
    reads = bam_file.fetch(good_chrom, begin, end, multiple_iterators=False)
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
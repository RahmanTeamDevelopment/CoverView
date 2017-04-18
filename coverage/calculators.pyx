"""
Utilities for calculating coverage summaries from a BAM file
"""

from __future__ import division


from coverage.statistics cimport QualityHistogramArray

from libc.stdint cimport uint32_t, uint64_t, uint8_t

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment, bam1_t, BAM_CIGAR_MASK,\
    BAM_CIGAR_SHIFT, BAM_CINS, BAM_CSOFT_CLIP, BAM_CREF_SKIP, BAM_CMATCH, BAM_CDEL, BAM_FDUP

import pysam

from pysam.libcalignmentfile cimport IteratorRowRegion
from cpython cimport array
from reads cimport ReadArray

import sys

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
    cdef array.array COV, QCOV, MEDBQ, FLBQ, MEDMQ, FLMQ
    cdef array.array COV_f, QCOV_f, MEDBQ_f, FLBQ_f, MEDMQ_f, FLMQ_f
    cdef array.array COV_r, QCOV_r, MEDBQ_r, FLBQ_r, MEDMQ_r, FLMQ_r
    cdef QualityHistogramArray bq_hists
    cdef QualityHistogramArray mq_hists

    def __init__(self, chrom, begin, end, bq_cutoff, mq_cutoff):
        cdef int bases_in_region = end - begin
        self.begin = begin
        self.end = end
        self.bq_cutoff = bq_cutoff
        self.mq_cutoff = mq_cutoff
        self.n_reads_in_region = 0
        self.COV = array.array('l', [0] * bases_in_region)
        self.QCOV = array.array('l', [0] * bases_in_region)
        self.MEDBQ = array.array('f', [float('NaN')] * bases_in_region)
        self.FLBQ = array.array('f', [float('NaN')] * bases_in_region)
        self.MEDMQ = array.array('f', [float('NaN')] * bases_in_region)
        self.FLMQ = array.array('f', [float('NaN')] * bases_in_region)

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
            self.FLMQ,
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


cdef void load_reads_into_array(ReadArray read_array, bam_file, chrom, start, end):
    """
    Load a chunk of BAM data into an in-memory read array
    """
    cdef int iterator_status = 0
    cdef IteratorRowRegion read_iterator = bam_file.fetch(chrom, start, end)

    while True:
        iterator_status = hts_itr_next(
            hts_get_bgzfp(read_iterator.htsfile),
            read_iterator.iter,
            read_iterator.b,
            read_iterator.htsfile
        )

        if iterator_status < 0:
            break

        read_array.append(read_iterator.b)


def get_profiles(bam_file, region, config):
    """
    Calculate and return coverage metrics for a specified region. Metrics include total
    coverage, coverage above the required base-quality and mapping quality threshold, fractions of
    low base qualities at a given position, fractions of low mapping quality reads covering a given
    position, median mapping qualities and median base qualities.

    This is by far the most computationally expensive part of CoverView. > 90% of the run-time is
    currently spent in this function.
    """
    cdef ReadArray read_array = ReadArray(100)

    bq_cutoff = float(config['low_bq'])
    mq_cutoff = float(config['low_mq'])

    chrom, beg_end = region.split(":")
    begin = int(beg_end.split("-")[0]) - 1
    end = int(beg_end.split("-")[1]) - 1 + 1 # TODO: sort out the indexing!

    good_chrom = get_valid_chromosome_name(chrom, bam_file)

    if not good_chrom in bam_file.references:
        return None

    load_reads_into_array(read_array, bam_file, good_chrom, begin, end)
    coverage_calc = SimpleCoverageCalculator(good_chrom, begin, end, bq_cutoff, mq_cutoff)
    reads = bam_file.fetch(good_chrom, begin, end, multiple_iterators=False)
    coverage_calc.add_reads(reads)

    return coverage_calc.get_coverage_summary()


def calculateChromData(samfile, ontarget):
    sys.stdout.write('\rFinalizing analysis ... 0.0%')
    sys.stdout.flush()

    chrom_lengths = samfile.lengths
    total_reference_length = sum(chrom_lengths)
    chroms = samfile.references

    chromdata = dict()
    chromsres = []
    total_mapped_reads_in_bam = 0
    allon = 0
    alloff = 0

    chrom_lengths_processed_so_far = 0

    for chrom,length in zip(chroms,chrom_lengths):

        total = get_num_mapped_reads_covering_chromosome(samfile, chrom)
        chrom_lengths_processed_so_far += length

        if 'chr' + chrom in ontarget:
            on = int(ontarget['chr' + chrom])
            off = total - on
        else:
            on = 0
            off = 0

        chromsres.append({'CHROM': chrom, 'RC': total, 'RCIN': on, 'RCOUT': off})
        total_mapped_reads_in_bam += total
        allon += on
        alloff += off

        x = round(100 * chrom_lengths_processed_so_far / total_reference_length, 1)
        sys.stdout.write('\rFinalizing analysis ... ' + str(x) + '%')
        sys.stdout.flush()

    chromdata['Chroms'] = chromsres
    chromdata['Mapped'] = {'RC': total_mapped_reads_in_bam, 'RCIN': allon, 'RCOUT': alloff}

    total_reads_in_bam = samfile.mapped + samfile.unmapped
    chromdata['Total'] = total_reads_in_bam
    chromdata['Unmapped'] = total_reads_in_bam - total_mapped_reads_in_bam

    sys.stdout.write('\rFinalizing analysis ... 100.0% - Done')
    sys.stdout.flush()
    print ''

    return chromdata


def calculateChromdata_minimal(samfile, options):
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

    allreads = pysam.flagstat(options.input)
    allreads = allreads[:allreads.find('+')]
    allreads = int(allreads.strip())

    chromdata['Total'] = allreads
    chromdata['Unmapped'] = allreads - alltotal

    # Finalizing progress info
    sys.stdout.write('\rRunning analysis ... 100.0% - Done')
    sys.stdout.flush()
    print ''

    return chromdata
"""
Utilities for calculating coverage summaries from a BAM file
"""

from __future__ import division

import tgmi.bamutils
import logging
import output

from cpython cimport array
from libc.stdint cimport uint32_t, uint64_t, uint8_t

from pysam.libcalignmentfile cimport IteratorRowRegion

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment, bam1_t, BAM_CIGAR_MASK,\
    BAM_CIGAR_SHIFT, BAM_CINS, BAM_CSOFT_CLIP, BAM_CREF_SKIP, BAM_CMATCH, BAM_CDEL, BAM_FDUP,\
    BAM_FREVERSE, bam_endpos

from pysam.libchtslib cimport BGZF, hts_get_bgzfp, hts_itr_t, hts_idx_t, htsFile, hts_itr_next, bam_get_qual,\
    bam_get_qname

from .reads cimport ReadArray
from .statistics cimport QualityHistogramArray

_logger = logging.getLogger("coverview")


cdef class RegionCoverageCalculator(object):
    """
    Utility class for computing coverage summaries for a specified genomic
    region
    """
    cdef int begin
    cdef int end
    cdef int bq_cutoff
    cdef int mq_cutoff
    cdef int n_reads_in_region, n_reads_in_region_f, n_reads_in_region_r
    cdef int count_duplicates
    cdef array.array COV, QCOV, MEDBQ, FLBQ, MEDMQ, FLMQ
    cdef array.array COV_f, QCOV_f, MEDBQ_f, FLBQ_f, MEDMQ_f, FLMQ_f
    cdef array.array COV_r, QCOV_r, MEDBQ_r, FLBQ_r, MEDMQ_r, FLMQ_r
    cdef QualityHistogramArray bq_hists, bq_hists_f, bq_hists_r
    cdef QualityHistogramArray mq_hists, mq_hists_f, mq_hists_r

    def __init__(self, chrom, begin, end, bq_cutoff, mq_cutoff, count_duplicates):
        cdef int bases_in_region = end - begin
        self.begin = begin
        self.end = end
        self.bq_cutoff = bq_cutoff
        self.mq_cutoff = mq_cutoff
        self.n_reads_in_region = 0
        self.n_reads_in_region_f = 0
        self.n_reads_in_region_r = 0

        self.COV = array.array('l', [0] * bases_in_region)
        self.COV_f = array.array('l', [0] * bases_in_region)
        self.COV_r = array.array('l', [0] * bases_in_region)

        self.QCOV = array.array('l', [0] * bases_in_region)
        self.QCOV_f = array.array('l', [0] * bases_in_region)
        self.QCOV_r = array.array('l', [0] * bases_in_region)

        self.MEDBQ = array.array('f', [float('NaN')] * bases_in_region)
        self.MEDBQ_f = array.array('f', [float('NaN')] * bases_in_region)
        self.MEDBQ_r = array.array('f', [float('NaN')] * bases_in_region)

        self.FLBQ = array.array('f', [float('NaN')] * bases_in_region)
        self.FLBQ_f = array.array('f', [float('NaN')] * bases_in_region)
        self.FLBQ_r = array.array('f', [float('NaN')] * bases_in_region)

        self.MEDMQ = array.array('f', [float('NaN')] * bases_in_region)
        self.MEDMQ_f = array.array('f', [float('NaN')] * bases_in_region)
        self.MEDMQ_r = array.array('f', [float('NaN')] * bases_in_region)

        self.FLMQ = array.array('f', [float('NaN')] * bases_in_region)
        self.FLMQ_f = array.array('f', [float('NaN')] * bases_in_region)
        self.FLMQ_r = array.array('f', [float('NaN')] * bases_in_region)

        self.bq_hists = QualityHistogramArray(bases_in_region)
        self.bq_hists_f = QualityHistogramArray(bases_in_region)
        self.bq_hists_r = QualityHistogramArray(bases_in_region)
        self.mq_hists = QualityHistogramArray(bases_in_region)
        self.mq_hists_f = QualityHistogramArray(bases_in_region)
        self.mq_hists_r = QualityHistogramArray(bases_in_region)

        if count_duplicates is True:
            self.count_duplicates = 1
        else:
            self.count_duplicates = 0

    cdef void add_reads(self, bam1_t** reads_start, bam1_t** reads_end):
        """
        """
        cdef int region_size = self.end - self.begin
        cdef int begin = self.begin
        cdef int end = self.end
        cdef int bq_cutoff = self.bq_cutoff
        cdef int mq_cutoff = self.mq_cutoff
        cdef long* COV = self.COV.data.as_longs
        cdef long* COV_f = self.COV_f.data.as_longs
        cdef long* COV_r = self.COV_r.data.as_longs
        cdef long* QCOV = self.QCOV.data.as_longs
        cdef long* QCOV_f = self.QCOV_f.data.as_longs
        cdef long* QCOV_r = self.QCOV_r.data.as_longs
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
        cdef int is_forward_read = 0

        while reads_start != reads_end:

            src = reads_start[0]

            if src.core.pos >= end:
                reads_start += 1
                continue

            if bam_endpos(src) <= begin:
                reads_start += 1
                continue

            # Skip duplicate reads
            if self.count_duplicates == 0 and src.core.flag & BAM_FDUP != 0:
                reads_start += 1
                continue

            # Reverse bit is set, so read is reverse
            if src.core.flag & BAM_FREVERSE != 0:
                is_forward_read = 0
            else:
                is_forward_read = 1

            n_cigar = src.core.n_cigar

            # This is an unmapped read whose mate is mapped in this region. Skip these.
            if n_cigar == 0:
                reads_start += 1
                continue

            # _logger.info("Read. Start = {}. End = {}. Name = {}. ISize = {}".format(
            #     src.core.pos,
            #     bam_endpos(src),
            #     bam_get_qname(src),
            #     src.core.isize
            # ))

            self.n_reads_in_region += 1

            if is_forward_read:
                self.n_reads_in_region_f += 1
            else:
                self.n_reads_in_region_r += 1

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
                        if begin <= i < end:
                            offset = i - begin
                            base_quality = base_qualities[index + (i-pos)]

                            self.bq_hists.add_data(offset, base_quality)
                            self.mq_hists.add_data(offset, mapping_quality)
                            COV[offset] += 1

                            if is_forward_read:
                                self.bq_hists_f.add_data(offset, base_quality)
                                self.mq_hists_f.add_data(offset, mapping_quality)
                                COV_f[offset] += 1
                            else:
                                self.bq_hists_r.add_data(offset, base_quality)
                                self.mq_hists_r.add_data(offset, mapping_quality)
                                COV_r[offset] += 1


                            if mapping_quality >= mq_cutoff and base_quality >= bq_cutoff:
                                QCOV[offset] += 1

                                if is_forward_read:
                                    QCOV_f[offset] += 1
                                else:
                                    QCOV_r[offset] += 1

                    pos += l
                    index += l

                elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                    for i from pos <= i < pos + l:
                        if begin <= i < end:
                            offset = i - begin
                            COV[offset] += 1

                            if is_forward_read:
                                COV_f[offset] += 1
                            else:
                                COV_r[offset] += 1

                            if mapping_quality >= mq_cutoff:
                                QCOV[offset] += 1

                                if is_forward_read:
                                    QCOV_f[offset] += 1
                                else:
                                    QCOV_r[offset] += 1

                    pos += l
            reads_start += 1

        self.compute_summary_statistics_for_region()

    cdef void compute_summary_statistics_for_region(self):
        cdef float* MEDBQ = self.MEDBQ.data.as_floats
        cdef float* MEDBQ_f = self.MEDBQ_f.data.as_floats
        cdef float* MEDBQ_r = self.MEDBQ_r.data.as_floats

        cdef float* MEDMQ = self.MEDMQ.data.as_floats
        cdef float* MEDMQ_f = self.MEDMQ_f.data.as_floats
        cdef float* MEDMQ_r = self.MEDMQ_r.data.as_floats

        cdef float* FLBQ = self.FLBQ.data.as_floats
        cdef float* FLBQ_f = self.FLBQ_f.data.as_floats
        cdef float* FLBQ_r = self.FLBQ_r.data.as_floats

        cdef float* FLMQ = self.FLMQ.data.as_floats
        cdef float* FLMQ_f = self.FLMQ_f.data.as_floats
        cdef float* FLMQ_r = self.FLMQ_r.data.as_floats

        for i in range(self.end - self.begin):
            MEDBQ[i] = self.bq_hists.compute_median(i)
            MEDBQ_f[i] = self.bq_hists_f.compute_median(i)
            MEDBQ_r[i] = self.bq_hists_r.compute_median(i)

            MEDMQ[i] = self.mq_hists.compute_median(i)
            MEDMQ_f[i] = self.mq_hists_f.compute_median(i)
            MEDMQ_r[i] = self.mq_hists_r.compute_median(i)

            FLBQ[i] = self.bq_hists.compute_fraction_below_threshold(i, <int>(self.bq_cutoff))
            FLBQ_f[i] = self.bq_hists_f.compute_fraction_below_threshold(i, <int>(self.bq_cutoff))
            FLBQ_r[i] = self.bq_hists_r.compute_fraction_below_threshold(i, <int>(self.bq_cutoff))

            FLMQ[i] = self.mq_hists.compute_fraction_below_threshold(i, <int>(self.mq_cutoff))
            FLMQ_f[i] = self.mq_hists_f.compute_fraction_below_threshold(i, <int>(self.mq_cutoff))
            FLMQ_r[i] = self.mq_hists_r.compute_fraction_below_threshold(i, <int>(self.mq_cutoff))

    def get_coverage_summary(self):
        return PerBaseCoverageSummary(
            self.n_reads_in_region,
            self.n_reads_in_region_f,
            self.n_reads_in_region_r,
            self.COV,
            self.QCOV,
            self.MEDBQ,
            self.FLBQ,
            self.MEDMQ,
            self.FLMQ,
            self.COV_f,
            self.QCOV_f,
            self.MEDBQ_f,
            self.FLBQ_f,
            self.MEDMQ_f,
            self.FLMQ_f,
            self.COV_r,
            self.QCOV_r,
            self.MEDBQ_r,
            self.FLBQ_r,
            self.MEDMQ_r,
            self.FLMQ_r
        )


cdef void load_reads_into_array(ReadArray read_array, bam_file, chrom, start, end):
    """
    Load a chunk of BAM data into an in-memory read array
    """
    cdef int iterator_status = 0

    _logger.info("Loading data for %s:%s-%s", chrom, start, end)

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


def get_region_coverage_summary(bam_file, cluster, config):
    """
    Calculate and return coverage metrics for a specified region. Metrics include total
    coverage, coverage above the required base-quality and mapping quality threshold, fractions of
    low base qualities at a given position, fractions of low mapping quality reads covering a given
    position, median mapping qualities and median base qualities.

    This is by far the most computationally expensive part of CoverView. > 90% of the run-time is
    currently spent in this function.
    
    """
    cdef ReadArray read_array = ReadArray(100)
    cdef bam1_t** reads_start
    cdef bam1_t** reads_end

    cluster_chrom = tgmi.bamutils.get_valid_chromosome_name(cluster[0].chromosome, bam_file)
    cluster_begin = cluster[0].start_pos
    cluster_end = cluster[-1].end_pos

    _logger.debug("Processing cluster of regions spanning {}:{}-{}".format(
        cluster_chrom, cluster_begin, cluster_end
    ))

    _logger.debug("Loading reads into in-memory array")

    load_reads_into_array(
        read_array,
        bam_file,
        cluster_chrom,
        cluster_begin,
        cluster_end
    )

    bq_cutoff = float(config['low_bq'])
    mq_cutoff = float(config['low_mq'])

    for interval in cluster:

        if config['outputs']['profiles'] or config['outputs']['regions']:

            coverage_calc = RegionCoverageCalculator(
                cluster_chrom,
                interval.start_pos,
                interval.end_pos,
                bq_cutoff,
                mq_cutoff,
                config['count_duplicate_reads']
            )

            read_array.set_pointers_to_start_and_end_of_interval(
                interval.start_pos,
                interval.end_pos,
                &reads_start,
                &reads_end
            )

            coverage_calc.add_reads(reads_start, reads_end)

            yield RegionCoverageSummary(
                interval.name,
                interval.chromosome,
                interval.start_pos,
                interval.end_pos,
                coverage_calc.get_coverage_summary()
            )


class BamFileCoverageSummary(object):
    """
    Overall coverage summary for one whole BAM file. Stores total numbers of reads,
    numbers of reads for each chromosome and other overall coverage summary metrics.
    """
    def __init__(
            self,
            num_reads_by_chromosome,
            num_reads,
            num_mapped_reads,
            num_on_target_reads,
            num_off_target_reads
    ):

        self.num_reads_by_chromosome = num_reads_by_chromosome
        self.num_reads = num_reads
        self.num_mapped_reads = num_mapped_reads
        self.num_on_target_reads = num_on_target_reads
        self.num_off_target_reads = num_off_target_reads

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """
        Return a string representation. This makes it possible to call print with
        instances of this class.
        """
        return str({
            "Chroms": self.num_reads_by_chromosome,
            "Mapped": {
                'RC': self.num_mapped_reads,
                'RCIN': self.num_on_target_reads,
                'RCOUT': self.num_off_target_reads
            },
            "Total": self.num_reads,
            "Unmapped": self.num_reads - self.num_mapped_reads
        })


class PerBaseCoverageSummary(object):
    """
    Stores summary coverage data for each base of a specified
    region, along with some overall region-level metrics.
    """
    def __init__(
            self,
            num_reads_in_region,
            num_forward_reads_in_region,
            num_reverse_reads_in_region,
            coverage_at_each_base,
            high_quality_coverage_at_each_base,
            median_quality_at_each_base,
            fraction_of_low_base_qualities_at_each_base,
            median_mapping_quality_at_each_base,
            fraction_of_low_mapping_qualities_at_each_base,
            forward_coverage_at_each_base,
            forward_high_quality_coverage_at_each_base,
            forward_median_quality_at_each_base,
            forward_fraction_of_low_base_qualities_at_each_base,
            forward_median_mapping_quality_at_each_base,
            forward_fraction_of_low_mapping_qualities_at_each_base,
            reverse_coverage_at_each_base,
            reverse_high_quality_coverage_at_each_base,
            reverse_median_quality_at_each_base,
            reverse_fraction_of_low_base_qualities_at_each_base,
            reverse_median_mapping_quality_at_each_base,
            reverse_fraction_of_low_mapping_qualities_at_each_base,
    ):
        self.num_reads_in_region = num_reads_in_region
        self.num_forward_reads_in_region = num_forward_reads_in_region
        self.num_reverse_reads_in_region = num_reverse_reads_in_region
        self.coverage_at_each_base = coverage_at_each_base
        self.high_quality_coverage_at_each_base = high_quality_coverage_at_each_base
        self.median_quality_at_each_base = median_quality_at_each_base
        self.fraction_of_low_base_qualities_at_each_base = fraction_of_low_base_qualities_at_each_base
        self.median_mapping_quality_at_each_base = median_mapping_quality_at_each_base
        self.fraction_of_low_mapping_qualities_at_each_base = fraction_of_low_mapping_qualities_at_each_base
        self.forward_coverage_at_each_base = forward_coverage_at_each_base
        self.forward_high_quality_coverage_at_each_base = forward_high_quality_coverage_at_each_base
        self.forward_median_quality_at_each_base = forward_median_quality_at_each_base
        self.forward_fraction_of_low_base_qualities_at_each_base = forward_fraction_of_low_base_qualities_at_each_base
        self.forward_median_mapping_quality_at_each_base = forward_median_mapping_quality_at_each_base
        self.forward_fraction_of_low_mapping_qualities_at_each_base = forward_fraction_of_low_mapping_qualities_at_each_base
        self.reverse_coverage_at_each_base = reverse_coverage_at_each_base
        self.reverse_high_quality_coverage_at_each_base = reverse_high_quality_coverage_at_each_base
        self.reverse_median_quality_at_each_base = reverse_median_quality_at_each_base
        self.reverse_fraction_of_low_base_qualities_at_each_base = reverse_fraction_of_low_base_qualities_at_each_base
        self.reverse_median_mapping_quality_at_each_base = reverse_median_mapping_quality_at_each_base
        self.reverse_fraction_of_low_mapping_qualities_at_each_base = reverse_fraction_of_low_mapping_qualities_at_each_base

    def __repr__(self):
        return self.__str__()

    def as_dict(self):
        return {
            "RC": self.num_reads_in_region,
            "RC_f": self.num_forward_reads_in_region,
            "RC_r": self.num_reverse_reads_in_region,
            "COV": list(self.coverage_at_each_base),
            "QCOV": list(self.high_quality_coverage_at_each_base),
            "MEDBQ": list(self.median_quality_at_each_base),
            "FLBQ":list( self.fraction_of_low_base_qualities_at_each_base),
            "MEDMQ": list(self.median_mapping_quality_at_each_base),
            "FLMQ": list(self.fraction_of_low_mapping_qualities_at_each_base),
            "COV_f": list(self.forward_coverage_at_each_base),
            "QCOV_f": list(self.forward_high_quality_coverage_at_each_base),
            "MEDBQ_f": list(self.forward_median_quality_at_each_base),
            "FLBQ_f": list(self.forward_fraction_of_low_base_qualities_at_each_base),
            "MEDMQ_f": list(self.forward_median_mapping_quality_at_each_base),
            "FLMQ_f": list(self.forward_fraction_of_low_mapping_qualities_at_each_base),
            "COV_r": list(self.reverse_coverage_at_each_base),
            "QCOV_r": list(self.reverse_high_quality_coverage_at_each_base),
            "MEDBQ_r": list(self.reverse_median_quality_at_each_base),
            "FLBQ_r": list(self.reverse_fraction_of_low_base_qualities_at_each_base),
            "MEDMQ_r": list(self.reverse_median_mapping_quality_at_each_base),
            "FLMQ_r": list(self.reverse_fraction_of_low_mapping_qualities_at_each_base)
        }

    def __str__(self):
        return str(self.as_dict())

    def print_to_file(
            self,
            bytes region_name,
            bytes chromosome,
            int start_position,
            int end_position,
            object transcript_database,
            int write_directional_summaries,
            object output_file,
            object low_quality_runs_output_file
    ):
        """
        Output the contents of this class to a text file.        
        """
        cdef array.array coverage_at_each_base = self.coverage_at_each_base
        cdef array.array high_quality_coverage_at_each_base = self.high_quality_coverage_at_each_base
        cdef array.array median_quality_at_each_base = self.median_quality_at_each_base
        cdef array.array fraction_of_low_base_qualities_at_each_base = self.fraction_of_low_base_qualities_at_each_base
        cdef array.array median_mapping_quality_at_each_base = self.median_mapping_quality_at_each_base
        cdef array.array fraction_of_low_mapping_qualities_at_each_base = self.fraction_of_low_mapping_qualities_at_each_base
        cdef array.array forward_coverage_at_each_base = self.forward_coverage_at_each_base
        cdef array.array forward_high_quality_coverage_at_each_base = self.forward_high_quality_coverage_at_each_base
        cdef array.array forward_median_quality_at_each_base = self.forward_median_quality_at_each_base
        cdef array.array forward_fraction_of_low_base_qualities_at_each_base = self.forward_fraction_of_low_base_qualities_at_each_base
        cdef array.array forward_median_mapping_quality_at_each_base = self.forward_median_mapping_quality_at_each_base
        cdef array.array forward_fraction_of_low_mapping_qualities_at_each_base = self.forward_fraction_of_low_mapping_qualities_at_each_base
        cdef array.array reverse_coverage_at_each_base = self.reverse_coverage_at_each_base
        cdef array.array reverse_high_quality_coverage_at_each_base = self.reverse_high_quality_coverage_at_each_base
        cdef array.array reverse_median_quality_at_each_base = self.reverse_median_quality_at_each_base
        cdef array.array reverse_fraction_of_low_base_qualities_at_each_base = self.reverse_fraction_of_low_base_qualities_at_each_base
        cdef array.array reverse_median_mapping_quality_at_each_base = self.reverse_median_mapping_quality_at_each_base
        cdef array.array reverse_fraction_of_low_mapping_qualities_at_each_base = self.reverse_fraction_of_low_mapping_qualities_at_each_base
        
        cdef long* COV_array = coverage_at_each_base.data.as_longs
        cdef long* QCOV_array = high_quality_coverage_at_each_base.data.as_longs
        cdef float* MEDBQ_array = median_quality_at_each_base.data.as_floats
        cdef float* FLBQ_array = fraction_of_low_base_qualities_at_each_base.data.as_floats
        cdef float* MEDMQ_array = median_mapping_quality_at_each_base.data.as_floats
        cdef float* FLMQ_array = fraction_of_low_mapping_qualities_at_each_base.data.as_floats
        cdef long* COV_f_array = forward_coverage_at_each_base.data.as_longs
        cdef long* QCOV_f_array = forward_high_quality_coverage_at_each_base.data.as_longs
        cdef float* MEDBQ_f_array = forward_median_quality_at_each_base.data.as_floats
        cdef float* FLBQ_f_array = forward_fraction_of_low_base_qualities_at_each_base.data.as_floats
        cdef float* MEDMQ_f_array = forward_median_mapping_quality_at_each_base.data.as_floats
        cdef float* FLMQ_f_array = forward_fraction_of_low_mapping_qualities_at_each_base.data.as_floats
        cdef long* COV_r_array = reverse_coverage_at_each_base.data.as_longs
        cdef long* QCOV_r_array = reverse_high_quality_coverage_at_each_base.data.as_longs
        cdef float* MEDBQ_r_array = reverse_median_quality_at_each_base.data.as_floats
        cdef float* FLBQ_r_array = reverse_fraction_of_low_base_qualities_at_each_base.data.as_floats
        cdef float* MEDMQ_r_array = reverse_median_mapping_quality_at_each_base.data.as_floats
        cdef float* FLMQ_r_array = reverse_fraction_of_low_mapping_qualities_at_each_base.data.as_floats

        cdef int num_bases = len(self.coverage_at_each_base)
        cdef int i, cov, qcov, cov_f, qcov_f, cov_r, qcov_r
        cdef float medbq, flbq, medmq, flmq
        cdef float medbq_f, flbq_f, medmq_f, flmq_f
        cdef float medbq_r, flbq_r, medmq_r, flmq_r
        cdef bytes output_line

        cdef int low_qual_window_start = -1
        cdef int low_qual_window_end = -1
        cdef bytes transcripts_overlapping_start_of_low_qual_window = None
        cdef bytes transcripts_overlapping_end_of_low_qual_window = None

        output_file.write('\n')
        output_file.write('[{}]\n'.format(region_name))

        for i from 0 <= i < num_bases:
            cov = COV_array[i]
            qcov = QCOV_array[i]
            medbq = MEDBQ_array[i]
            flbq = FLBQ_array[i]
            medmq = MEDMQ_array[i]
            flmq = FLMQ_array[i]
            cov_f = COV_f_array[i]
            qcov_f = QCOV_f_array[i]
            medbq_f = MEDBQ_f_array[i]
            flbq_f = FLBQ_f_array[i]
            medmq_f = MEDMQ_f_array[i]
            flmq_f = FLMQ_f_array[i]
            cov_r = COV_r_array[i]
            qcov_r = QCOV_r_array[i]
            medbq_r = MEDBQ_r_array[i]
            flbq_r = FLBQ_r_array[i]
            medmq_r = MEDMQ_r_array[i]
            flmq_r = FLMQ_r_array[i]

            if write_directional_summaries == 0:
                output_line = "{}\t{}\t{}\t{}\t{}\t{:.3f}\t{}\t{:.3f}\n".format(
                    chromosome,
                    start_position + i,
                    cov,
                    qcov,
                    medbq,
                    flbq,
                    medmq,
                    flmq
                )

                output_line = output_line.replace("nan", ".")
                output_file.write(output_line)
            else:
                output_line = \
                    "{}\t{}\t{}\t{}\t{}\t{:.3f}\t{}\t{:.3f}\t{}\t{}\t{}\t{:.3f}\t{}\t{:.3f}" \
                    "\t{}\t{}\t{}\t{:.3f}\t{}\t{:.3f}\n" \
                    "".format(
                        chromosome,
                        start_position + i,
                        cov,
                        qcov,
                        medbq,
                        flbq,
                        medmq,
                        flmq,
                        cov_f,
                        qcov_f,
                        medbq_f,
                        flbq_f,
                        medmq_f,
                        flmq_f,
                        cov_r,
                        qcov_r,
                        medbq_r,
                        flbq_r,
                        medmq_r,
                        flmq_r
                    )

                output_line = output_line.replace("nan", ".")
                output_file.write(output_line)

            if transcript_database is not None:
                if qcov < 15:
                    if low_qual_window_start == -1:
                        low_qual_window_start = start_position + i
                        transcripts_overlapping_start_of_low_qual_window = output.get_transcripts_overlapping_position(
                            transcript_database,
                            chromosome,
                            start_position + i
                        )
                else:
                    if low_qual_window_start != -1:
                        transcripts_overlapping_end_of_low_qual_window = output.get_transcripts_overlapping_position(
                            transcript_database,
                            chromosome,
                            start_position + i - 1
                        )

                    low_quality_runs_output_file.write(
                        "{}\t{]\t{}\t{}\t{}\t{}\n".format(
                            region_name,
                            chromosome,
                            low_qual_window_start,
                            low_qual_window_end,
                            transcripts_overlapping_start_of_low_qual_window,
                            transcripts_overlapping_end_of_low_qual_window
                        )
                    )


class RegionCoverageSummary(object):
    """
    Stores data summarising coverage for a genomic region
    """
    def __init__(
            self,
            region_name,
            chromosome,
            start_position,
            end_position,
            per_base_coverage_profile
    ):
        self.region_name = region_name
        self.chromosome = chromosome
        self.start_position = start_position
        self.end_position = end_position
        self.per_base_coverage_profile = per_base_coverage_profile

    def as_dict(self):
        return {
            "Name": self.region_name,
            "Chrom": self.chromosome,
            "Start": self.start_position,
            "End": self.end_position,
            "profiles": self.per_base_coverage_profile.as_dict()
        }

    def __str__(self):
        return str(self.as_dict)


def calculate_chromosome_coverage_metrics(bam_file, on_target):
    _logger.info("Calculating per-chromosome coverage metrics")

    chromosomes = bam_file.references
    chromosome_lengths = bam_file.lengths

    number_of_reads_covering_chromosomes = []

    total_reads_in_bam = bam_file.mapped + bam_file.unmapped
    total_on_target_reads = 0
    total_off_target_reads = 0

    bam_index_stats = tgmi.bamutils.load_bam_index_stats_from_file(bam_file)

    for chrom, length in zip(chromosomes, chromosome_lengths):

        num_mapped_reads = bam_index_stats.get_num_mapped_reads_for_chromosome(chrom)

        if chrom in on_target:
            num_on_target_reads = on_target[chrom]
            num_off_target_reads = num_mapped_reads - num_on_target_reads
        else:
            num_on_target_reads = 0
            num_off_target_reads = num_mapped_reads

        number_of_reads_covering_chromosomes.append({
            'CHROM': chrom,
            'RC': num_mapped_reads,
            'RCIN': num_on_target_reads,
            'RCOUT': num_off_target_reads
        })

        total_on_target_reads += num_on_target_reads
        total_off_target_reads += num_off_target_reads

    _logger.info("Finished calculating per-chromosome coverage metrics")

    return {
        "Chroms": number_of_reads_covering_chromosomes,
        "Mapped": {
            'RC': bam_index_stats.get_total_mapped_reads_in_bam(),
            'RCIN': total_on_target_reads,
            'RCOUT': total_off_target_reads
        },
        "Total": bam_index_stats.get_total_reads_in_bam(),
        "Unmapped": bam_index_stats.get_total_unmapped_reads_in_bam()
    }


def calculate_minimal_chromosome_coverage_metrics(bam_file, options):
    _logger.info("Calculating minimal per-chromosome coverage metrics")

    bam_index_stats = tgmi.bamutils.load_bam_index_stats_from_file(bam_file)
    number_of_reads_covering_chromosomes = []

    for chrom in bam_file.references:
        num_reads = bam_index_stats.get_num_total_reads_for_chromosome(chrom)
        number_of_reads_covering_chromosomes.append({
            'CHROM': chrom, 'RC': num_reads
        })

    _logger.info("Finished calculating minimal per-chromosome coverage metrics")

    return {
        "Chroms": number_of_reads_covering_chromosomes,
        "Mapped": {
            "RC": bam_index_stats.get_total_mapped_reads_in_bam()
        },
        "Total": bam_index_stats.get_total_reads_in_bam(),
        "Unmapped": bam_index_stats.get_total_unmapped_reads_in_bam()
    }

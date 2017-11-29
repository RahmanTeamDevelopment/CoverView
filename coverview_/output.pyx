"""
Format and write per-base coverage profile output
"""

import csv
import logging

from . import transcript


_logger = logging.getLogger("coverview_")
_canonical_chromosomes = set( range(1, 23) + ["X", "Y", "MT"] )


def get_transcripts_overlapping_position(transcript_database, chrom, pos):
    """
    Returns a comma-separated list of the transcripts which overlap this base, and the coordinate of
    the base within each transcript, in.
    """
    transcript_coordinates = transcript.get_transcript_coordinates(transcript_database, chrom, pos)
    transcripts = []

    for key, value in transcript_coordinates.items():
        transcripts.append(key.gene_symbol + ':' + key.ensembl_id + ':' + value)

    if len(transcripts) == 0:
        return '.'

    return ','.join(transcripts)


class PerBaseCoverageOutput(object):
    """
    Data and functions needed for producing summary coverage output for
    each base across a region.
    """
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.out_profiles = open(options.output + '_profiles.txt', 'w')
        self.out_poor = None
        self.only_output_profiles_for_flagged_regions = False
        self.output_directional_coverage_summaries = config['direction']
        self.output_directional_coverage_information = False

        if config['only_flagged_profiles']:
            self.only_output_profiles_for_flagged_regions = True

        if options.transcript_db is not None and config['transcript'].get('poor') is True:
            self.out_poor = open(options.output + '_poor.txt', 'w')

        if config['direction']:
            self.output_directional_coverage_information = True

    def __del__(self):
        self.out_profiles.close()

        if self.out_poor is not None:
            self.out_poor.close()

    def write_header(self):

        profheader = [
            'Chromosome',
            'Position'
        ]

        if self.options.transcript_db is not None and self.config['transcript'].get('profiles') is True:
            _logger.info("Transcript coordinates will be written in profiles output")
            profheader.append('Transcript_coordinate')

        profheader.extend([
            'COV',
            'QCOV',
            'MEDBQ',
            'FLBQ',
            'MEDMQ',
            'FLMQ'
        ])

        if self.output_directional_coverage_information:
            profheader.extend([
                'COV+',
                'QCOV+',
                'MEDBQ+',
                'FLBQ+',
                'MEDMQ+',
                'FLMQ+',
                'COV-',
                'QCOV-',
                'MEDBQ-',
                'FLBQ-',
                'MEDMQ-',
                'FLMQ-'
            ])

        self.out_profiles.write('#' + '\t'.join(profheader) + '\n')

        if self.out_poor is not None:
            poorheader = [
                'Region',
                'Chromosome',
                'Start_position',
                'End_position',
                'Start_transcript',
                'End_transcript'
            ]

            self.out_poor.write('#' + '\t'.join(poorheader) + '\n')

    def write_output(self, coverage_data, transcript_database):
        if self.only_output_profiles_for_flagged_regions and coverage_data.passes_thresholds:
            pass
        else:
            per_base_summary = coverage_data.per_base_coverage_profile

            if self.options.transcript_db is not None and self.config['transcript'].get('profiles') is True:
                write_transcripts_in_profiles = 1
            else:
                write_transcripts_in_profiles = 0

            per_base_summary.print_to_file(
                coverage_data.region_name,
                coverage_data.chromosome,
                coverage_data.start_position,
                coverage_data.end_position,
                transcript_database,
                self.output_directional_coverage_information,
                self.out_profiles,
                self.out_poor,
                write_transcripts_in_profiles
            )


class RegionsOutput(object):
    """
    Data and functions needed for producing summary coverage output for
    targeted regions. There are several optional metrics, which will only be output
    if the relevant configuration options are set.
    """
    def __init__(self, options, config):
        self.output_file = open(options.output + '_regions.txt', 'w')
        self.config = config
        self.options = options

        if options.transcript_db is not None and config['transcript'].get('regions') is True:
            self.output_transcript_data = True
        else:
            self.output_transcript_data = False

        if not self.config['pass'] is None:
            self.output_region_pass_flag = True
        else:
            self.output_region_pass_flag = False

        if self.config['direction']:
            self.output_directional_coverage_information = True
        else:
            self.output_directional_coverage_information = False

    def __del__(self):
        self.output_file.close()

    def write_header(self):
        """
        Write the file header. This is done only once. The header has several
        optional fields (e.g. for directional coverage summaries) which will only
        be output if the relevant configuration parameters were set.
        """
        header = [
            'Region',
            'Chromosome',
            'Start_position',
            'End_position'
        ]

        if self.output_transcript_data:
            header.extend([
                'Start_transcript',
                'End_transcript'
            ])

        if self.output_region_pass_flag:
            header.append(
                'Pass_or_flag'
            )

        header.extend([
            'RC',
            'MEDCOV',
            'MINCOV',
            'MEDQCOV',
            'MINQCOV',
            'MAXFLMQ',
            'MAXFLBQ'
        ])

        if self.output_directional_coverage_information:
            header.extend([
                'RC+',
                'MEDCOV+',
                'MINCOV+',
                'MEDQCOV+',
                'MINQCOV+',
                'MAXFLMQ+',
                'MAXFLBQ+',
                'RC-',
                'MEDCOV-',
                'MINCOV-',
                'MEDQCOV-',
                'MINQCOV-',
                'MAXFLMQ-',
                'MAXFLBQ-'
            ])

        self.output_file.write(
            '#' + '\t'.join(header) + '\n'
        )

    def write_output(self, coverage_data, transcript_database):

        region_name = coverage_data.region_name
        chrom = coverage_data.chromosome
        region_start = coverage_data.start_position
        region_end = coverage_data.end_position
        coverage_summary = coverage_data.summary

        if not chrom.startswith("chr") and chrom in _canonical_chromosomes:
            chrom = "chr{}".format(chrom)

        output_record = [
            region_name,
            chrom,
            region_start,
            region_end
        ]

        if self.output_transcript_data:

            transcripts_overlapping_start = get_transcripts_overlapping_position(
                transcript_database, chrom, region_start
            )

            transcripts_overlapping_end = get_transcripts_overlapping_position(
                transcript_database, chrom, region_end
            )

            output_record.extend([
                transcripts_overlapping_start,
                transcripts_overlapping_end
            ])

        if self.output_region_pass_flag:
            if coverage_data.passes_thresholds:
                output_record.append('PASS')
            else:
                output_record.append('FLAG')

        output_record.extend([
            coverage_data.per_base_coverage_profile.num_reads_in_region,
            coverage_summary['MEDCOV'],
            coverage_summary['MINCOV'],
            coverage_summary['MEDQCOV'],
            coverage_summary['MINQCOV'],
            coverage_summary['MAXFLMQ'],
            coverage_summary['MAXFLBQ']
        ])

        if self.output_directional_coverage_information:
            output_record.extend([
                coverage_data.per_base_coverage_profile.num_forward_reads_in_region,
                coverage_summary['MEDCOV_f'],
                coverage_summary['MINCOV_f'],
                coverage_summary['MEDQCOV_f'],
                coverage_summary['MINQCOV_f'],
                coverage_summary['MAXFLMQ_f'],
                coverage_summary['MAXFLBQ_f'],
                coverage_data.per_base_coverage_profile.num_reverse_reads_in_region,
                coverage_summary['MEDCOV_r'],
                coverage_summary['MINCOV_r'],
                coverage_summary['MEDQCOV_r'],
                coverage_summary['MINQCOV_r'],
                coverage_summary['MAXFLMQ_r'],
                coverage_summary['MAXFLBQ_r']
            ])

        self.output_file.write(
            ('\t'.join(str(x) for x in output_record) + '\n').replace("nan", ".")
        )


def output_chromosome_coverage_metrics(options, chromosome_coverage_metrics):
    """
    Write the per-chromosome coverage metrics to a tab-separated file    
    """
    num_mapped_reads_in_bam = chromosome_coverage_metrics['Mapped']
    num_unmapped_reads_in_bam = chromosome_coverage_metrics['Unmapped']
    num_total_reads_in_bam = chromosome_coverage_metrics['Total']

    with open(options.output + '_summary.txt', 'wb') as output_file:
        csv_writer = csv.writer(output_file, delimiter='\t')

        csv_writer.writerows([
            ["#CHROM", "RC", "RCIN", "RCOUT"],
            ["Total", num_total_reads_in_bam, "-", "-"],
            ["Unmapped", num_unmapped_reads_in_bam, "-", "-"],
            ["Mapped", num_mapped_reads_in_bam["RC"], num_mapped_reads_in_bam["RCIN"], num_mapped_reads_in_bam["RCOUT"]]
        ])

        for coverage_metrics in chromosome_coverage_metrics['Chroms']:
            chromosome = coverage_metrics['CHROM']

            csv_writer.writerow([
                chromosome,
                coverage_metrics["RC"],
                coverage_metrics["RCIN"],
                coverage_metrics["RCOUT"]
            ])


def output_minimal_chromosome_coverage_metrics(options, chromosome_coverage_metrics):
    """
    Write the minimal (i.e. in/out counts for targeted regions) per-chromosome
    coverage metrics to a tab-separated file    
    """
    num_mapped_reads_in_bam = chromosome_coverage_metrics['Mapped']
    num_unmapped_reads_in_bam = chromosome_coverage_metrics['Unmapped']
    num_total_reads_in_bam = chromosome_coverage_metrics['Total']

    with open(options.output + '_summary.txt', 'wb') as output_file:
        csv_writer = csv.writer(output_file, delimiter='\t')

        csv_writer.writerows([
            ["#CHROM", "RC"],
            ["Total", num_total_reads_in_bam],
            ["Unmapped", num_unmapped_reads_in_bam],
            ["Mapped", num_mapped_reads_in_bam["RC"]]
        ])

        for coverage_metrics in chromosome_coverage_metrics['Chroms']:
            chromosome = coverage_metrics['CHROM']
            csv_writer.writerow([
                chromosome,
                coverage_metrics["RC"]
            ])


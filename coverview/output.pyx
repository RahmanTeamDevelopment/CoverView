"""
Format and write per-base coverage profile output
"""

import csv
import datetime
import json
import logging
import pysam
import transcript


_logger = logging.getLogger("coverview")
_canonical_chromosomes = set( range(1, 23) + ["X", "Y", "MT"] )


def get_transcripts_overlapping_position(transcript_database, chrom, pos):
    """
    Returns a comma-separated list of the transcripts which overlap this base, and the coordinate of
    the base within each transcript, in.
    """
    transcript_coordinates = transcript.getTranscriptCoordinates(transcript_database, chrom, pos)
    transcripts = []

    for key, value in transcript_coordinates.iteritems():
        transcripts.append(key.geneSymbol + ':' + key.ENST + ':' + value)

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
        self.only_output_profiles_for_failed_regions = False
        self.output_directional_coverage_summaries = config['direction']
        self.output_directional_coverage_information = False

        if config['only_fail_profiles']:
            self.only_output_profiles_for_failed_regions = True

        if config['transcript_db'] is not None:
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

        if self.out_poor is not None:
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

    def write_output(self, coverage_data):
        if self.only_output_profiles_for_failed_regions and coverage_data.passes_thresholds:
            pass
        else:
            per_base_summary = coverage_data.per_base_coverage_profile

            per_base_summary.print_to_file(
                coverage_data.region_name,
                coverage_data.chromosome,
                coverage_data.start_position,
                coverage_data.end_position,
                self.config['transcript_db'],
                self.output_directional_coverage_information,
                self.out_profiles,
                self.out_poor
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

        if config['transcript_db'] is not None:
            self.output_transcript_data = True
        else:
            self.output_transcript_data = False

        if not self.config['pass'] is None:
            self.output_region_pass_fail = True
        else:
            self.output_region_pass_fail = False

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

        if self.output_region_pass_fail:
            header.append(
                'Pass_or_fail'
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
                'MEDCOV+',
                'MINCOV+',
                'MEDQCOV+',
                'MINQCOV+',
                'MAXFLMQ+',
                'MAXFLBQ+',
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

    def write_output(self, coverage_data):

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
                self.enstdb, chrom, region_start
            )

            transcripts_overlapping_end = get_transcripts_overlapping_position(
                self.entsdb, chrom, region_end
            )

            output_record.extend([
                transcripts_overlapping_start,
                transcripts_overlapping_end
            ])

        if self.output_region_pass_fail:
            if coverage_data.passes_thresholds:
                output_record.append('PASS')
            else:
                output_record.append('FAIL')

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
                coverage_summary['MEDCOV_f'],
                coverage_summary['MINCOV_f'],
                coverage_summary['MEDQCOV_f'],
                coverage_summary['MINQCOV_f'],
                coverage_summary['MAXFLMQ_f'],
                coverage_summary['MAXFLBQ_f'],
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


class GuiOutput(object):
    """
    Encapsulates all data and functions needed or producing the JSON output for
    the graphical interface.
    """
    def __init__(self, options, config):
        self.output_file = open(options.output  + '_gui/data/results.js', 'w')
        self.ref_file = pysam.Fastafile(config["reference"])
        self.have_written_first_line = False

    def __del__(self):
        self.output_file.write(']')
        self.output_file.close()
        self.ref_file.close()

    def write_header(self):
        self.output_file.write('function readData() {\n')
        self.output_file.write('\tdata={\"targets\":[')

    def write_output(self, coverage_data):

        self.output_json(coverage_data)

        if not self.have_written_first_line:
            self.have_written_first_line = True

    def output_json(self, coverage_data):
        if self.have_written_first_line:
            self.out_json.write(',')

        self.output_file.write(
            json.dumps(coverage_data, separators=(',', ':'))
        )

    def finalize_output(self, options, chromdata, config,
                        numOfTargets, failedtargets, uniqueIDs, uniqueids):

        newchromsres = []

        others = {
            'CHROM': '...',
            'RC': 0,
            'RCIN': 0,
            'RCOUT': 0
        }

        chrnames = [
            '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19',
            '20', '21', '22', 'X', 'Y', 'chr1', 'chr2', 'chr3', 'chr4',
            'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
            'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
            'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
        ]

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

        self.output_file.write(
            ',\"chromdata\":' + json.dumps(chromdata, separators=(',', ':'))
        )

        infn = options.input

        if '/' in infn:
            infn = infn[infn.rfind('/') + 1:]

        if len(infn) > 27:
            infn = infn[:27] + '...bam'

        self.output_file.write(',\"input\":\"' + infn + '\"')
        self.output_file.write(',\"direction\":' + str(config['direction']).lower())
        self.output_file.write(',\"duplicates\":' + str(config['duplicates']).lower())
        self.output_file.write(',\"ntargets\":' + str(numOfTargets))
        self.output_file.write(',\"unique\":' + str(uniqueIDs))
        self.output_file.write(',\"nfailed\":' + str(failedtargets))
        self.output_file.write(',\"uniquefailed\":' + str(len(uniqueids)))

        bedfn = options.bedfile

        if '/' in bedfn:
            bedfn = bedfn[bedfn.rfind('/') + 1:]

        if len(bedfn) > 20:
            bedfn = bedfn[:20] + '...bed'

        self.output_file.write(',\"bedfile\":\"' + bedfn + '\"')

        now = datetime.datetime.now()
        self.output_file.write(',\"date\":\"' + now.strftime("%d-%m-%Y, %H:%M") + "\"")

        passdef = []

        for k, v in config['pass'].iteritems():
            [met,minmax] = k.split('_')
            if minmax == 'MIN':
                passdef.append(met + '>' + str(v))
            else:
                passdef.append(met + '<' + str(v))

        passdefstr = ', '.join(passdef)
        self.output_file.write(',\"passdef\":\"' + passdefstr + '\"')

        passmets = dict()

        for k, v in config['pass'].iteritems():
            [met, minmax] = k.split('_')
            if met.endswith('QCOV'):
                passmets['QCOV'] = v
            elif met.endswith('COV'):
                passmets['COV'] = v
            elif met.endswith('FLMQ'):
                passmets['FLMQ'] = v
            elif met.endswith('FLBQ'):
                passmets['FLBQ'] = v
            passmets[met] = v

        self.output_file.write(',\"passmets\":' + json.dumps(passmets, separators=(',', ':')))
        self.output_file.write('}\n')
        self.output_file.write('\treturn data\n')
        self.output_file.write('}\n')


def output_chromosome_coverage_metrics(options, chromosome_coverage_metrics):
    """
    Write the per-chromosome coverage metrics to a tab-separated file    
    """
    num_mapped_reads_in_bam = chromosome_coverage_metrics['Mapped']
    num_unmapped_reads_in_bam = chromosome_coverage_metrics['Unmapped']
    num_total_reads_in_bam = chromosome_coverage_metrics['Total']

    _logger.info(num_mapped_reads_in_bam)
    _logger.info(num_unmapped_reads_in_bam)
    _logger.info(num_total_reads_in_bam)

    with open(options.output + '_summary.txt', 'w') as output_file:
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

    with open(options.output + '_summary.txt', 'w') as output_file:
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


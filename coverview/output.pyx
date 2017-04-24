"""
Format and write per-base coverage profile output
"""

import datetime
import json
import logging

from cpython cimport array

import cython
import pysam
cimport cython


logger = logging.getLogger("coverview")


class ProfilesOutput(object):
    """
    Data and functions needed for producing summary coverage output for
    each base across a region.
    """
    def __init__(self, options, config):
        pass

    def __del__(self):
        pass

    def write_header(self):
        pass

    def compute_coverage_metric(self, coverage_data):
        pass


class RegionsOutput(object):
    """
    Data and functions needed for producing summary coverage output for
    targeted regions.
    """
    def __init__(self, options, config):
        self.output_file = open(options.output + '_regions.txt', 'w')
        self.config = config
        self.options = options

    def __del__(self):
        self.out_targets.close()

    def write_header(self):
        header = ['Region', 'Chromosome', 'Start_position', 'End_position']

        if self.config['transcript_db'] is not None:
            header.extend(
                ['Start_transcript', 'End_transcript']
            )

        if self.config['pass'] is not None:
            header.append('Pass_or_fail')

        header.extend(
            ['RC', 'MEDCOV', 'MINCOV', 'MEDQCOV', 'MINQCOV', 'MAXFLMQ', 'MAXFLBQ']
        )

        if self.config['direction']:
            header.extend(
                ['MEDCOV+', 'MINCOV+', 'MEDQCOV+', 'MINQCOV+', 'MAXFLMQ+', 'MAXFLBQ+']
            )
            header.extend(
                ['MEDCOV-', 'MINCOV-', 'MEDQCOV-', 'MINQCOV-', 'MAXFLMQ-', 'MAXFLBQ-']
            )

        self.output_file.write('#' + '\t'.join(header) + '\n')

    def compute_coverage_metric(self, coverage_data):
        coverage_data['Summary'] = self.compute_summaries_of_region_coverage(
            coverage_data['Profiles']
        )

        if self.config['pass'] is not None:
            coverage_data['PASS'] = self.targetPASS(coverage_data)
        else:
            coverage_data['PASS'] = True

        if self.config['outputs']['gui']:
            coverage_data['Ref'] = self.getReferenceSequence(chrom, begin, end)
        else:
            coverage_data['Ref'] = ''


class GuiOutput(object):
    """
    Encapsulates all data and functions needed or producing the JSON output for
    the graphical interface.
    """
    def __init__(self, options, config):
        self.output_file = open(options.output  + '_gui/data/results.js', 'w')
        self.ref_file = pysam.Fastafile(config["reference"])

    def __del__(self):
        self.output_file.write(']')
        self.output_file.close()
        self.ref_file.close()

    def write_header(self):
        self.output_file.write('function readData() {\n')
        self.output_file.write('\tdata={\"targets\":[')

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


def output_target_file_header(config, out_poor, out_json, out_targets, out_profiles):

    if config['outputs']['profiles']:
        profheader = ['Chromosome', 'Position']

        if config['transcript']['profiles'] and not config['transcript_db'] is None:
            profheader.append('Transcript_coordinate')

        profheader.extend(
            ['COV', 'QCOV', 'MEDBQ', 'FLBQ', 'MEDMQ', 'FLMQ']
        )

        if config['direction']:
            profheader.extend(
                ['COV+', 'QCOV+', 'MEDBQ+', 'FLBQ+', 'MEDMQ+', 'FLMQ+']
            )
            profheader.extend(
                ['COV-', 'QCOV-', 'MEDBQ-', 'FLBQ-', 'MEDMQ-', 'FLMQ-']
            )

        out_profiles.write('#' + '\t'.join(profheader) + '\n')

        if not config['transcript_db'] is None and config['outputs']['profiles']:
            poorheader = [
                'Region', 'Chromosome', 'Start_position',
                'End_position', 'Start_transcript',
                'End_transcript'
            ]

            out_poor.write('#' + '\t'.join(poorheader) + '\n')

        if config['outputs']['gui']:



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


def output_summary_minimal(options, chromdata):
    out_summary = open(
        options.output + '_summary.txt', 'w'
    )

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


# def output_profiles_with_transcript_coordinates(config, target, enstdb, output_file):
#     profiles = target['Profiles']
#
#     output_file.write('\n')
#     output_file.write('[{}]\n'.format(target['Name']))
#
#     if config['transcript_db'] is not None:
#         window_qcov = {
#             "targetname": None,
#             "chrom": None,
#             "start": None,
#             "transcriptstart": None
#         }
#
#     num_bases = len(profiles['COV'])
#
#     COV = profiles['COV']
#     QCOV = profiles['QCOV']
#     MEDBQ = profiles['MEDBQ']
#     FLBQ = profiles['FLBQ']
#     MEDMQ = profiles['MEDMQ']
#     FLMQ = profiles['FLMQ']
#
#     chrom = target['Chrom']
#     target_start_pos = target['Start']
#
#     for i in xrange(num_bases):
#         transcoordstr = ''
#
#         if config['transcript']['profiles'] and not config['transcript_db'] is None:
#
#             transcoords = transcript.getTranscriptCoordinates(
#                 enstdb,
#                 chrom,
#                 target_start_pos + i
#             )
#
#             v = []
#
#             for key, value in transcoords.iteritems():
#                 v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
#
#             transcoordstr = ','.join(v)
#
#         if config['transcript']['profiles'] and not config['transcript_db'] is None:
#             record.append(transcoordstr)
#
#         if config['transcript_db'] is not None:
#             if profiles['QCOV'][i] < 15:
#                 if window_qcov['transcriptstart'] is None:
#                     window_qcov['targetname'] = target['Name']
#                     window_qcov['chrom'] = target['Chrom']
#                     window_qcov['start'] = target['Start'] + i
#
#                     if transcoordstr == '':
#                         transcoords = transcript.getTranscriptCoordinates(
#                             enstdb,
#                             target['Chrom'],
#                             target['Start'] + i
#                         )
#
#                         v = []
#                         for key, value in transcoords.iteritems():
#                             v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
#
#                         transcoordstr = ','.join(v)
#                     window_qcov['transcriptstart'] = transcoordstr
#             else:
#                 if not window_qcov['transcriptstart'] is None:
#
#                     transcoords = transcript.getTranscriptCoordinates(
#                         enstdb,
#                         target['Chrom'],
#                         target['Start'] + i - 1
#                     )
#                     v = []
#
#                     for key, value in transcoords.iteritems():
#                         v.append(
#                             key.geneSymbol + ':' + key.ENST + ':' + value
#                         )
#
#                     transcoordstr_end = ','.join(v)
#
#                     record = [
#                         window_qcov['targetname'],
#                         window_qcov['chrom'],
#                         str(window_qcov['start']),
#                         str(target['Start'] + i - 1),
#                         window_qcov['transcriptstart'],
#                         transcoordstr_end
#                     ]
#
#                     if record[4] == '':
#                         record[4] = 'None'
#
#                     if record[5] == '':
#                         record[5] = 'None'
#
#                     if not (record[4] == 'None' and record[5] == 'None'):
#                         out_poor.write('\t'.join(record) + '\n')
#
#                     window_qcov = {
#                         "targetname": None,
#                         "chrom": None,
#                         "start": None,
#                         "transcriptstart": None
#                     }
#
#     if config['transcript_db'] is not None:
#         if not window_qcov['transcriptstart'] is None:
#
#             if transcoordstr == '':
#                 transcoords = transcript.getTranscriptCoordinates(
#                     enstdb,
#                     target['Chrom'],
#                     target['Start'] + i
#                 )
#
#                 v = []
#
#                 for key, value in transcoords.iteritems():
#                     v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
#
#                 transcoordstr = ','.join(v)
#
#             record = [
#                 window_qcov['targetname'],
#                 window_qcov['chrom'],
#                 str(window_qcov['start']),
#                 str(target['Start'] + i), window_qcov['transcriptstart'],
#                 transcoordstr
#             ]
#
#             if record[4] == '':
#                 record[4] = 'None'
#
#             if record[5] == '':
#                 record[5] = 'None'
#
#             if not record[4] == 'None' and record[5] == 'None':
#                 output_file.write('\t'.join(record) + '\n')


@cython.profile(False)
cdef bytes format_output_line(bytes chrom,
                              int start_pos,
                              int cov,
                              int qcov,
                              float medbq,
                              float flbq,
                              float medmq,
                              float flmq):
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            chrom,
            start_pos,
            cov,
            qcov,
            medbq,
            flbq,
            medmq,
            flmq
        )


def output_profiles(dict target, output_file):
    cdef dict profiles = target['Profiles']
    cdef array.array COV = profiles['COV']
    cdef array.array QCOV = profiles['QCOV']
    cdef array.array MEDBQ = profiles['MEDBQ']
    cdef array.array FLBQ = profiles['FLBQ']
    cdef array.array MEDMQ = profiles['MEDMQ']
    cdef array.array FLMQ = profiles['FLMQ']

    cdef long* COV_array = COV.data.as_longs
    cdef long* QCOV_array = QCOV.data.as_longs
    cdef float* MEDBQ_array = MEDBQ.data.as_floats
    cdef float* FLBQ_array = FLBQ.data.as_floats
    cdef float* MEDMQ_array = MEDMQ.data.as_floats
    cdef float* FLMQ_array = FLMQ.data.as_floats

    output_file.write('\n')
    output_file.write('[{}]\n'.format(target['Name']))

    cdef bytes chrom = target['Chrom']
    cdef int target_start_pos = target['Start']

    cdef int num_bases = len(COV)
    cdef int i, cov, qcov
    cdef float medbq, flbq, medmq, flmq
    cdef bytes output_line

    for i from 0 <= i < num_bases:
        cov = COV_array[i]
        qcov = QCOV_array[i]
        medbq = MEDBQ_array[i]
        flbq = FLBQ_array[i]
        medmq = MEDMQ_array[i]
        flmq = FLMQ_array[i]

        output_line = format_output_line(
            chrom,
            target_start_pos + i,
            cov,
            qcov,
            medbq,
            flbq,
            medmq,
            flmq)

        output_line = output_line.replace("nan", ".")
        output_file.write(output_line)


# def output_directional_profiles(config, target, enstdb, output_file):
#     profiles = target['Profiles']
#
#     output_file.write('\n')
#     output_file.write('[{}]\n'.format(target['Name']))
#
#     if config['transcript_db'] is not None:
#         window_qcov = {
#             "targetname": None,
#             "chrom": None,
#             "start": None,
#             "transcriptstart": None
#         }
#
#     num_bases = len(profiles['COV'])
#
#     COV = profiles['COV']
#     QCOV = profiles['QCOV']
#     MEDBQ = profiles['MEDBQ']
#     FLBQ = profiles['FLBQ']
#     MEDMQ = profiles['MEDMQ']
#     FLMQ = profiles['FLMQ']
#
#     chrom = target['Chrom']
#     target_start_pos = target['Start']
#
#     for i in xrange(num_bases):
#         transcoordstr = ''
#
#         if config['transcript']['profiles'] and not config['transcript_db'] is None:
#
#             transcoords = transcript.getTranscriptCoordinates(
#                 enstdb,
#                 chrom,
#                 target_start_pos + i
#             )
#
#             v = []
#
#             for key, value in transcoords.iteritems():
#                 v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
#
#             transcoordstr = ','.join(v)
#
#         record = [
#             chrom,
#             str(target_start_pos + i)
#         ]
#
#         if config['transcript']['profiles'] and not config['transcript_db'] is None:
#             record.append(transcoordstr)
#
#         record.extend(
#             [
#                 COV[i],
#                 QCOV[i],
#                 MEDBQ[i],
#                 FLBQ[i],
#                 MEDMQ[i],
#                 FLMQ[i]
#             ]
#         )
#
#         if config['direction']:
#             record.extend(
#                 [
#                     str(profiles['COV_f'][i]),
#                     str(profiles['QCOV_f'][i]),
#                     str(profiles['MEDBQ_f'][i]),
#                     str(profiles['FLBQ_f'][i]),
#                     str(profiles['MEDMQ_f'][i]),
#                     str(profiles['FLMQ_f'][i]),
#                     str(profiles['COV_r'][i]),
#                     str(profiles['QCOV_r'][i]),
#                     str(profiles['MEDBQ_r'][i]),
#                     str(profiles['FLBQ_r'][i]),
#                     str(profiles['MEDMQ_r'][i]),
#                     str(profiles['FLMQ_r'][i])
#                 ]
#             )
#
#         if record[-1] == 'nan':
#             record[-1] = '.'
#
#         if record[-2] == 'nan':
#             record[-2] = '.'
#
#         if record[-3] == 'nan':
#             record[-3] = '.'
#
#         if record[-4] == 'nan':
#             record[-4] = '.'
#
#         output_file.write('\t'.join([str(x) for x in record]) + '\n')
#
#         if config['transcript_db'] is not None:
#             if profiles['QCOV'][i] < 15:
#                 if window_qcov['transcriptstart'] is None:
#                     window_qcov['targetname'] = target['Name']
#                     window_qcov['chrom'] = target['Chrom']
#                     window_qcov['start'] = target['Start'] + i
#
#                     if transcoordstr == '':
#                         transcoords = transcript.getTranscriptCoordinates(
#                             enstdb,
#                             target['Chrom'],
#                             target['Start'] + i
#                         )
#
#                         v = []
#                         for key, value in transcoords.iteritems():
#                             v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
#
#                         transcoordstr = ','.join(v)
#                     window_qcov['transcriptstart'] = transcoordstr
#             else:
#                 if not window_qcov['transcriptstart'] is None:
#
#                     transcoords = transcript.getTranscriptCoordinates(
#                         enstdb,
#                         target['Chrom'],
#                         target['Start'] + i - 1
#                     )
#                     v = []
#
#                     for key, value in transcoords.iteritems():
#                         v.append(
#                             key.geneSymbol + ':' + key.ENST + ':' + value
#                         )
#
#                     transcoordstr_end = ','.join(v)
#
#                     record = [
#                         window_qcov['targetname'],
#                         window_qcov['chrom'],
#                         str(window_qcov['start']),
#                         str(target['Start'] + i - 1),
#                         window_qcov['transcriptstart'],
#                         transcoordstr_end
#                     ]
#
#                     if record[4] == '':
#                         record[4] = 'None'
#
#                     if record[5] == '':
#                         record[5] = 'None'
#
#                     if not (record[4] == 'None' and record[5] == 'None'):
#                         out_poor.write('\t'.join(record) + '\n')
#
#                     window_qcov = {
#                         "targetname": None,
#                         "chrom": None,
#                         "start": None,
#                         "transcriptstart": None
#                     }
#
#     if config['transcript_db'] is not None:
#         if not window_qcov['transcriptstart'] is None:
#
#             if transcoordstr == '':
#                 transcoords = transcript.getTranscriptCoordinates(
#                     enstdb,
#                     target['Chrom'],
#                     target['Start'] + i
#                 )
#
#                 v = []
#
#                 for key, value in transcoords.iteritems():
#                     v.append(key.geneSymbol + ':' + key.ENST + ':' + value)
#
#                 transcoordstr = ','.join(v)
#
#             record = [
#                 window_qcov['targetname'],
#                 window_qcov['chrom'],
#                 str(window_qcov['start']),
#                 str(target['Start'] + i), window_qcov['transcriptstart'],
#                 transcoordstr
#             ]
#
#             if record[4] == '':
#                 record[4] = 'None'
#
#             if record[5] == '':
#                 record[5] = 'None'
#
#             if not record[4] == 'None' and record[5] == 'None':
#                 output_file.write('\t'.join(record) + '\n')
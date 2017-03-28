"""
Format and write per-base coverage profile output
"""

from cpython cimport array


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
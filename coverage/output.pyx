"""
Format and write per-base coverage profile output
"""

from cpython cimport array


def output_target_file_header(self, config, out_poor, out_json, out_targets, out_profiles):
    if config['outputs']['regions']:

        targetheader = ['Region', 'Chromosome', 'Start_position', 'End_position']

        if config['transcript']['regions'] and not config['transcript_db'] is None:
            targetheader.extend(
                ['Start_transcript', 'End_transcript']
            )

        if not config['pass'] is None:
            targetheader.append('Pass_or_fail')

        targetheader.extend(
            ['RC', 'MEDCOV', 'MINCOV', 'MEDQCOV', 'MINQCOV', 'MAXFLMQ', 'MAXFLBQ']
        )

        if config['direction']:
            targetheader.extend(
                ['MEDCOV+', 'MINCOV+', 'MEDQCOV+', 'MINQCOV+', 'MAXFLMQ+', 'MAXFLBQ+']
            )
            targetheader.extend(
                ['MEDCOV-', 'MINCOV-', 'MEDQCOV-', 'MINQCOV-', 'MAXFLMQ-', 'MAXFLBQ-']
            )

        out_targets.write('#' + '\t'.join(targetheader) + '\n')

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
            out_json.write('function readData() {\n')
            out_json.write('\tdata={\"targets\":[')


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

    if config['transcript_db'] is not None:
        if  (config['outputs']['regions'] and config['transcript']['regions']) or (config['outputs']['profiles'] and config['transcript']['profiles']):
            print "Transcript db file:     " + config['transcript_db']
            print ''

    formats = 'summary'
    if config['outputs']['regions']: formats += ', regions'
    if config['outputs']['profiles']:
        if config['only_fail_profiles']:
            formats += ', profiles (failed regions)'
        else:
            formats += ', profiles (all regions)'
    if config['transcript_db'] is not None and config['outputs']['profiles']:
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
    print "--------------------------------------------------------------------------------------"


def printInfo_minimal(options):
    print 'Input, output and settings:'
    print "--------------------------------------------------------------------------------------"
    print "Input file name:        " + options.input
    print ''
    print "Output formats:         _summary"
    print "Output files prefix:    " + options.output
    print "--------------------------------------------------------------------------------------"


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


def finalizeJSONOutput(options):

    out_json = open(options.output + '_gui/data/results.js', 'a')

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

    out_json.write(
        ',\"chromdata\":' + json.dumps(chromdata, separators=(',', ':'))
    )

    infn = options.input

    if '/' in infn:
        infn = infn[infn.rfind('/') + 1:]

    if len(infn) > 27:
        infn = infn[:27] + '...bam'

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

    if len(bedfn) > 20:
        bedfn = bedfn[:20] + '...bed'

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
        if met.endswith('QCOV'):
            passmets['QCOV'] = v
        elif met.endswith('COV'):
            passmets['COV'] = v
        elif met.endswith('FLMQ'):
            passmets['FLMQ'] = v
        elif met.endswith('FLBQ'):
            passmets['FLBQ'] = v
        passmets[met] = v

    out_json.write(',\"passmets\":' + json.dumps(passmets, separators=(',', ':')))
    out_json.write('}\n')
    out_json.write('\treturn data\n')
    out_json.write('}\n')
    out_json.close()


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
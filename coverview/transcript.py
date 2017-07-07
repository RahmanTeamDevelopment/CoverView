from collections import OrderedDict


class Transcript(object):
    def __init__(self, line):
        self.exons = []
        cols = line.split('\t')
        self.ensembl_id = cols[0]
        self.gene_symbol = cols[1]
        self.gene_ID = cols[2]
        self.chrom = cols[4]
        self.strand = int(cols[5])
        self.transcript_start = int(cols[6])
        self.transcript_end = int(cols[7])
        self.coding_start = int(cols[8])
        self.coding_start_genomic = int(cols[9])
        self.coding_end_genomic = int(cols[10])

        for i in range(1, len(cols) - 11, 2):
            self.exons.append(
                Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i]))
            )

    def is_in_utr(self, pos):
        if self.strand == 1:
            return (pos < self.coding_start_genomic) or (pos > self.coding_end_genomic)
        else:
            return (pos > self.coding_start_genomic) or (pos < self.coding_end_genomic)


class Exon(object):
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start


def get_transcript_coordinates(enstdb, chrom, pos):
    ret = OrderedDict()
    transcripts = find_transcripts(enstdb, chrom, pos)

    for ensembl_id, transcript in transcripts.items():
        x, y = transform_to_csn_coordinate(pos, transcript)
        transcript_coordinates = 'c.' + str(x)
        if y != 0:
            if y > 0:
                transcript_coordinates += '+' + str(y)
            else:
                transcript_coordinates += str(y)
        ret[transcript] = transcript_coordinates
    return ret


def find_transcripts(enstdb, chrom, pos):
    ret = OrderedDict()

    if chrom in enstdb.contigs:
        good_chrom = chrom
    else:
        if 'chr' + chrom in enstdb.contigs:
            good_chrom = 'chr' + chrom
        else:
            if chrom.startswith('chr') and chrom[3:] in enstdb.contigs:
                good_chrom = chrom[3:]
            else:
                return ret

    # Checking both end points of the variant
    reg = good_chrom + ':' + str(pos) + '-' + str(pos)
    hits = enstdb.fetch(region=reg)

    for line in hits:
        transcript = Transcript(line)
        if not transcript.transcript_start + 1 <= pos <= transcript.transcript_end:
            continue
        ret[transcript.ensembl_id] = transcript

    return ret


def transform_to_csn_coordinate(pos, transcript):
    previous_exon_end = 99999999

    if not transcript.is_in_utr(pos):
        sum_of_exon_lengths = -transcript.codingStart + 1

        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]
            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if previous_exon_end < pos < exon.start + 1:
                        if pos <= (exon.start + 1 - previous_exon_end) / 2 + previous_exon_end:
                            x, y = transform_to_csn_coordinate(previous_exon_end, transcript)
                            return x, pos - previous_exon_end
                        else:
                            x, y = transform_to_csn_coordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < previous_exon_end:
                        if pos >= (previous_exon_end - exon.end + 1) / 2 + exon.end:
                            x, y = transform_to_csn_coordinate(previous_exon_end, transcript)
                            return x, previous_exon_end - pos
                        else:
                            x, y = transform_to_csn_coordinate(exon.end, transcript)
                            return x, exon.end - pos

            # Checking if genomic position is within exon
            if exon.start + 1 <= pos <= exon.end:
                if transcript.strand == 1:
                    return sum_of_exon_lengths + pos - exon.start, 0
                else:
                    return sum_of_exon_lengths + exon.end - pos + 1, 0

            # Calculating sum of exon lengths up to this point
            sum_of_exon_lengths += exon.length
            if transcript.strand == 1:
                previous_exon_end = exon.end
            else:
                previous_exon_end = exon.start + 1

    # If genomic position is within UTR
    else:
        if transcript.strand == 1:
            if pos < transcript.coding_start_genomic:
                return pos - transcript.coding_start_genomic, 0
            if pos > transcript.coding_end_genomic:
                return '+' + str(pos - transcript.coding_end_genomic), 0
        else:
            if pos > transcript.coding_start_genomic:
                return transcript.coding_start_genomic - pos, 0
            if pos < transcript.coding_end_genomic:
                return '+' + str(transcript.coding_end_genomic - pos), 0

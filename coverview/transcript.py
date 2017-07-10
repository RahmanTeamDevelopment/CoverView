import pysam
from collections import OrderedDict


CHROM_IDX = 2
TRANSCRIPT_START_IDX = 6
TRANSCRIPT_END_IDX = 7


def write_transcripts_to_indexed_tabix_file(transcripts, file_name):
    with open(file_name, 'w') as transcript_database:
        for transcript in transcripts:
            transcript_database.write(
                convert_transcript_to_line_of_ensemble_database(
                    transcript
                )
            )

    pysam.tabix_compress(
        file_name,
        file_name + '.gz',
        force=True
    )

    pysam.tabix_index(
        file_name + '.gz',
        seq_col=CHROM_IDX,
        start_col=TRANSCRIPT_START_IDX,
        end_col=TRANSCRIPT_END_IDX,
        meta_char='#',
        force=True
    )


def create_transcript_from_line_of_ensembl_database(line):
    cols = line.split('\t')
    exons = []
    ensembl_id = cols[0]
    gene_symbol = cols[1]
    gene_id = cols[2]
    chrom = cols[4]
    strand = int(cols[5])
    transcript_start = int(cols[6])
    transcript_end = int(cols[7])
    coding_start = int(cols[8])
    coding_start_genomic = int(cols[9])
    coding_end_genomic = int(cols[10])

    for i in range(1, len(cols) - 11, 2):
        exons.append(
            Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i]))
        )

    return Transcript(
        ensembl_id=ensembl_id,
        gene_symbol=gene_symbol,
        gene_id=gene_id,
        chrom=chrom,
        strand=strand,
        transcript_start=transcript_start,
        transcript_end=transcript_end,
        coding_start=coding_start,
        coding_start_genomic=coding_start_genomic,
        coding_end_genomic=coding_end_genomic,
        exons=exons
    )


def convert_transcript_to_line_of_ensemble_database(transcript):
    cols = [
        transcript.ensembl_id,
        transcript.gene_symbol,
        transcript.gene_id,
        "DUMMY_VALUE",
        transcript.chrom,
        str(transcript.strand),
        str(transcript.transcript_start),
        str(transcript.transcript_end),
        str(transcript.coding_start),
        str(transcript.coding_start_genomic),
        str(transcript.coding_end_genomic),
    ]

    for exon in transcript.exons:
        cols.append(str(exon.start))
        cols.append(str(exon.end))

    return "\t".join(cols)


class Transcript(object):
    def __init__(
            self,
            ensembl_id,
            gene_symbol,
            gene_id,
            chrom,
            strand,
            transcript_start,
            transcript_end,
            coding_start,
            coding_start_genomic,
            coding_end_genomic,
            exons
    ):
        self.ensembl_id = ensembl_id
        self.gene_symbol = gene_symbol
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.transcript_start = transcript_start
        self.transcript_end = transcript_end
        self.coding_start = coding_start
        self.coding_start_genomic = coding_start_genomic
        self.coding_end_genomic = coding_end_genomic
        self.exons = exons

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


def get_transcript_coordinates(transcript_database, chrom, pos):
    ret = OrderedDict()
    transcripts = find_transcripts(transcript_database, chrom, pos)

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


def find_transcripts(transcript_database, chrom, pos):
    ret = OrderedDict()

    if chrom in transcript_database.contigs:
        good_chrom = chrom
    else:
        if 'chr' + chrom in transcript_database.contigs:
            good_chrom = 'chr' + chrom
        else:
            if chrom.startswith('chr') and chrom[3:] in transcript_database.contigs:
                good_chrom = chrom[3:]
            else:
                return ret

    # Checking both end points of the variant
    reg = good_chrom + ':' + str(pos) + '-' + str(pos)
    hits = transcript_database.fetch(region=reg)

    for line in hits:

        transcript = create_transcript_from_line_of_ensembl_database(
            line
        )

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

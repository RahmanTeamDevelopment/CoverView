import logging
import pysam

from collections import OrderedDict


_logger = logging.getLogger("coverview")


CHROM_IDX = 4
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
        assert transcript_start == exons[0].start
        assert transcript_end == exons[-1].end

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
            return (pos < self.coding_start_genomic) or (pos >= self.coding_end_genomic)
        else:
            return (pos > self.coding_start_genomic) or (pos < self.coding_end_genomic)


class Exon(object):
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start

    def __contains__(self, position):
        return self.start <= position < self.end


def find_transcripts(transcript_database, chrom, pos):
    overlapping_transcripts = OrderedDict()
    region = "{}:{}-{}".format(chrom, pos, pos+1)
    hits = transcript_database.fetch(region=region)

    for line in hits:
        transcript = create_transcript_from_line_of_ensembl_database(
            line
        )

        if transcript.transcript_start <= pos < transcript.transcript_end:
            overlapping_transcripts[transcript.ensembl_id] = transcript

    return overlapping_transcripts


def get_transcript_coordinates(transcript_database, chrom, pos):
    ret = OrderedDict()
    transcripts = find_transcripts(transcript_database, chrom, pos)

    for ensembl_id, transcript in transcripts.items():
        transcript_coordinates = get_csn_coordinates(
            pos,
            transcript
        )

        ret[transcript] = transcript_coordinates
    return ret


def get_csn_coordinates(position, transcript):
    x, y = transform_to_csn_coordinate(position, transcript)
    transcript_coordinates = 'c.' + str(x)

    if y != 0:
        if y > 0:
            transcript_coordinates += '+' + str(y)
        else:
            transcript_coordinates += str(y)

    return transcript_coordinates


def transform_to_csn_coordinate(pos, transcript):

    if transcript.is_in_utr(pos):
        if pos < transcript.coding_start_genomic:
            return "{}".format(transcript.strand * (pos - transcript.coding_start_genomic)), 0
        elif pos >= transcript.coding_end_genomic:
            return "+{}".format(transcript.strand * (pos - (transcript.coding_end_genomic -1))), 0
        else:
            raise StandardError("This should not happen.")

    previous_exon = None
    sum_of_exon_lengths = -transcript.coding_start

    for i in range(len(transcript.exons)):
        exon = transcript.exons[i]
        if i > 0:
            if transcript.strand == 1:
                # Checking if genomic position is within intron
                if previous_exon.end <= pos < exon.start:
                    if pos <= (exon.start + 1 - previous_exon.end) / 2 + previous_exon.end:
                        return sum_of_exon_lengths, pos - (previous_exon.end - 1)
                    else:
                        return sum_of_exon_lengths + 1, pos - exon.start
            else:
                # Checking if genomic position is within intron
                if exon.end <= pos < previous_exon.start:
                    if pos >= (previous_exon.start - exon.end + 1) / 2 + exon.end:
                        return sum_of_exon_lengths + 1, previous_exon.start - pos
                    else:
                        return sum_of_exon_lengths, exon.end - pos

        # Checking if genomic position is within exon
        if pos in exon:
            if transcript.strand == 1:
                return sum_of_exon_lengths + pos - exon.start + 1, 0
            else:
                return sum_of_exon_lengths + exon.end - pos + 1, 0

        sum_of_exon_lengths += exon.length
        previous_exon = exon



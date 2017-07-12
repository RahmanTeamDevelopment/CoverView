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

        if strand == 1:
            assert transcript_start == exons[0].start
            assert transcript_end == exons[-1].end
        else:
            assert transcript_start == exons[-1].start
            assert transcript_end == exons[0].end

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

    def is_position_in_utr(self, pos):
        if self.strand == 1:
            return (pos < self.coding_start_genomic) or (pos >= self.coding_end_genomic)
        else:
            return (pos > self.coding_start_genomic) or (pos <= self.coding_end_genomic)

    def get_distance_from_coding_region(self, position):
        if self.strand == 1:
            if position < self.coding_start_genomic:
                return position - self.coding_start_genomic
            elif position >= self.coding_end_genomic:
                return position - (self.coding_end_genomic - 1)
            else:
                return 0
        else:
            if position > self.coding_start_genomic:
                return self.coding_start_genomic - position
            elif position <= self.coding_end_genomic:
                return (self.coding_end_genomic - 1) - position
            else:
                return 0


class Exon(object):
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start

    def __contains__(self, position):
        return self.start <= position < self.end

    def distance_from(self, position):
        if position < self.start:
            return position - self.start
        elif position >= self.end:
            return position - (self.end - 1)
        else:
            return 0

    def get_coordinate_of_position_forward(self, position):
        """
        Returns the 0-based offset in this exon of the position when traversing the exon in the
        forward direction.
        """
        if position not in self:
            return None
        else:
            return position - self.start

    def get_coordinate_of_position_reverse(self, position):
        """
        Returns the 0-based offset in this exon of the position when traversing the exon in the
        reverse direction.
        """
        if position not in self:
            return None
        else:
            return self.end - position - 1


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

    if transcript.is_position_in_utr(position):
        distance_from_coding_region = transcript.get_distance_from_coding_region(position)
        return "c.{0:+d}".format(distance_from_coding_region)
    else:
        coding_position, distance_to_exon = get_position_in_coding_sequence(
            position,
            transcript
        )

        transcript_coordinates = 'c.{}'.format(coding_position)

        if distance_to_exon != 0:
            transcript_coordinates += '{0:+d}'.format(distance_to_exon)

        return transcript_coordinates


def get_position_in_coding_sequence(position, transcript):

    # We always report 1-based positions in CSN
    position_in_coding_sequence = 1 - transcript.coding_start
    previous_exon = None

    for exon in transcript.exons:
        if previous_exon is not None:
            distance_from_left_exon = previous_exon.distance_from(position)
            distance_from_right_exon = exon.distance_from(position)

            if transcript.strand == 1:

                if distance_from_left_exon != 0 and distance_from_right_exon != 0:
                    if abs(distance_from_left_exon) <= abs(distance_from_right_exon):
                        return position_in_coding_sequence - 1, distance_from_left_exon
                    else:
                        return position_in_coding_sequence, distance_from_right_exon
            else:
                if distance_from_left_exon != 0 and distance_from_right_exon != 0:
                    if abs(distance_from_right_exon) <= abs(distance_from_left_exon):
                        return position_in_coding_sequence, distance_from_right_exon
                    else:
                        return position_in_coding_sequence - 1, distance_from_left_exon

        if position in exon:
            if transcript.strand == 1:
                position_in_exon = exon.get_coordinate_of_position_forward(position)
            else:
                position_in_exon = exon.get_coordinate_of_position_reverse(position)

            return position_in_coding_sequence + position_in_exon, 0

        # This takes us to the position of the first base in the next exon, so we need to -1 when referring back
        # to the previous exon for intronic coordinates.
        position_in_coding_sequence += exon.length
        previous_exon = exon


def old_transform_to_csn_coordinate(pos, transcript):
    prev_exon_end = 99999999
    # Checking if genomic position is within translated region
    if not transcript.is_position_in_utr(pos):
        sum_of_exon_lengths = -transcript.coding_start + 1
        # Iterating through exons
        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]
            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prev_exon_end < pos < exon.start + 1:
                        if pos <= (exon.start + 1 - prev_exon_end) / 2 + prev_exon_end:
                            x, y = old_transform_to_csn_coordinate(prev_exon_end, transcript)
                            return x, pos - prev_exon_end
                        else:
                            x, y = old_transform_to_csn_coordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prev_exon_end:
                        if pos >= (prev_exon_end - exon.end + 1) / 2 + exon.end:
                            x, y = old_transform_to_csn_coordinate(prev_exon_end, transcript)
                            return x, prev_exon_end - pos
                        else:
                            x, y = old_transform_to_csn_coordinate(exon.end, transcript)
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
                prev_exon_end = exon.end
            else:
                prev_exon_end = exon.start + 1
    # If genomic position is within UTR
    else:
        if transcript.strand == 1:
            if pos < transcript.coding_start_genomic: return pos - transcript.coding_start_genomic, 0
            if pos > transcript.coding_end_genomic: return '+' + str(pos - transcript.coding_end_genomic), 0
        else:
            if pos > transcript.coding_start_genomic: return transcript.coding_start_genomic - pos, 0
            if pos < transcript.coding_end_genomic: return '+' + str(transcript.coding_end_genomic - pos), 0

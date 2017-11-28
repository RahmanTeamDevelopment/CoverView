import logging
import pysam

from collections import OrderedDict


_logger = logging.getLogger("coverview_")


CHROM_IDX = 4
TRANSCRIPT_START_IDX = 6
TRANSCRIPT_END_IDX = 7


def write_transcripts_to_indexed_tabix_file(transcripts, file_name):
    with open(file_name, 'w') as transcript_database:
        for transcript in transcripts:
            transcript_database.write(
                convert_transcript_to_line_of_old_database(
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


def create_transcript_from_line_of_old_database(line):
    """
    The old database format used a mixture of 0-based and 1-based coordinates, with some half-open intervals
    some closed. This is being replaced with a newer format that uses 0-based half-open intervals everywhere.
    The assertions here check the data according to the old format, in order to prevent this function being used
    with the wrong data. Once everything is ported over to the new data, this function should be updated or replaced.
    """
    cols = line.split('\t')
    exons = []
    ensembl_id = cols[0]
    gene_symbol = cols[1]
    gene_id = cols[2]
    chrom = cols[4]
    strand = int(cols[5])
    transcript_start = int(cols[6])
    transcript_end = int(cols[7])
    coding_start = int(cols[8]) - 1
    coding_start_genomic = int(cols[9]) - 1

    for i in range(1, len(cols) - 11, 2):
        exons.append(
            Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i]))
        )

    if strand == 1:
        # The value for 1-based inclusive is the same as for 0-based half-open
        coding_end_genomic = int(cols[10])
    else:
        # To convert to 0-based, half-open in reverse direction
        coding_end_genomic = int(cols[10]) - 2

    # If the first exon is short then these assertion will not hold as the coding_start may be in the next
    # exon with an intron between
    if strand == 1:
        if exons[0].length >= coding_start and not transcript_start + coding_start == coding_start_genomic:
            _logger.error("Invalid forward transcript data in input database")
            _logger.error(line)
    else:
        if exons[0].length >= coding_start and not transcript_end - 1 - coding_start == coding_start_genomic:
            _logger.error("Invalid reverse transcript data in input database")
            _logger.error(line)

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


def convert_transcript_to_line_of_old_database(transcript):
    """
    As with the above function, this should be removed / updated when we move to the new database
    format.
    """
    if transcript.strand == 1:
        cols = [
            transcript.ensembl_id,
            transcript.gene_symbol,
            transcript.gene_id,
            "DUMMY_VALUE",
            transcript.chrom,
            str(transcript.strand),
            str(transcript.transcript_start),
            str(transcript.transcript_end),
            str(transcript.coding_start + 1),
            str(transcript.coding_start_genomic + 1),
            str(transcript.coding_end_genomic),
        ]
    else:
        cols = [
            transcript.ensembl_id,
            transcript.gene_symbol,
            transcript.gene_id,
            "DUMMY_VALUE",
            transcript.chrom,
            str(transcript.strand),
            str(transcript.transcript_start),
            str(transcript.transcript_end),
            str(transcript.coding_start + 1),
            str(transcript.coding_start_genomic + 1),
            str(transcript.coding_end_genomic + 2)
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

        one_based_position, distance_to_exon = get_position_in_coding_sequence(
            self.coding_end_genomic,
            self
        )

        self.total_length_of_coding_sequence = one_based_position - 1


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

    try:
        hits = transcript_database.fetch(region=region)
    except ValueError:
        return overlapping_transcripts

    for line in hits:
        transcript = create_transcript_from_line_of_old_database(
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

    coding_position, distance_to_exon = get_position_in_coding_sequence(
        position,
        transcript
    )

    if coding_position <= 0:
        coding_position -= 1

    if coding_position > transcript.total_length_of_coding_sequence:
        transcript_coordinates = 'c.+{}'.format(
            coding_position - transcript.total_length_of_coding_sequence
        )
    else:
        transcript_coordinates = 'c.{}'.format(coding_position)

    if distance_to_exon != 0:
        transcript_coordinates += '{0:+d}'.format(distance_to_exon)

    return transcript_coordinates


def get_position_in_coding_sequence(position, transcript):

    # We always report 1-based positions in CSN, hence the 1 here rather than just
    # - transcript.coding start which would be correct for 0-indexing.
    position_in_coding_sequence = 1 - transcript.coding_start
    previous_exon = None

    for exon in transcript.exons:
        if previous_exon is not None:
            distance_from_previous_exon = previous_exon.distance_from(position)
            distance_from_exon = exon.distance_from(position)

            if transcript.strand == 1:

                if distance_from_previous_exon > 0 > distance_from_exon:
                    if abs(distance_from_previous_exon) <= abs(distance_from_exon):
                        return position_in_coding_sequence - 1, distance_from_previous_exon
                    else:
                        return position_in_coding_sequence, distance_from_exon
            else:
                if distance_from_previous_exon < 0 < distance_from_exon:
                    if abs(distance_from_exon) <= abs(distance_from_previous_exon):
                        return position_in_coding_sequence, -1 * distance_from_exon
                    else:
                        return position_in_coding_sequence - 1, -1 * distance_from_previous_exon

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

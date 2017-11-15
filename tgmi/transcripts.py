"""Module for working with transcripts, creating transcript DBs and reading from transcript DB"""

from __future__ import division
import gzip
from operator import itemgetter
import pysam
import datetime
import os


class Transcript(object):
    """Class for a single transcript"""


    def __init__(self, id=None, version=None, hgnc_id=None, chrom=None, start=None, end=None, strand=None, exons=None,
                 coding_start=None, coding_end=None, cdna_length=None, prot_length=None, info=None):
        """Constructor of Transcript class

            Notes:
                    start, end: 0-based coordinates coordinates, inclusive for start, exclusive for end
                    coding_start, coding_end: 0-based coordinates (inclusive for both)
                    strand: + or -
                    exons: list of Exon objects
        """

        self.id = id
        self.version = version
        self.hgnc_id = hgnc_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = exons
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.cdna_length = cdna_length
        self.prot_length = prot_length
        self.info = info

    def read_from_database_record(self, record):
        """Read transcript from transcript database record

            Notes:
                record is a dict of column:value pairs
        """

        # Set attributes of transcript
        for key, value in record.iteritems():
            if key in ['start', 'end', 'coding_start', 'coding_end', 'cdna_coding_start', 'cdna_coding_end']:
                value = int(value)
            if key == 'exons':
                self.exons = [Exon(s) for s in value.split(',')]
                continue
            setattr(self, key, value)

        # Finalize transcript
        self.finalize()

    def set_info(self):
        """Calculate and sets the info field"""

        if self._any_unset(['strand', 'start', 'end', 'exons', 'cdna_length', 'prot_length']): return
        self.info = '+/' if self.strand == '+' else '-/'
        self.info += str(round((self.end - self.start + 1) / 1000, 1)) + 'kb/'
        self.info += str(len(self.exons)) + '/' + str(round(self.cdna_length / 1000, 1)) + 'kb/' + str(self.prot_length)

    def cds_regions(self):
        """Return list of CDS regions

            Notes:
                Returned is a list of (x,y) pairs where coordinates x and y are 0-based, inclusive for x and exclusive for y
        """

        if self._any_unset(['exons', 'coding_start', 'coding_end']):
            return

        ret = []
        for exon in self.exons:
            cds_region = exon.get_cds(self.coding_start, self.coding_end)
            if cds_region is not None:
                ret.append(cds_region)
        return ret

    def utr5_regions(self):
        """Return list of UTR5 regions

            Notes:
                Returned is a list of (x,y) pairs where coordinates x and y are 0-based, inclusive for x and exclusive for y
        """

        if self._any_unset(['exons', 'coding_start', 'strand']):
            return

        ret = []
        for exon in self.exons:
            utr5_region = exon.get_utr5(self.coding_start, self.strand)
            if utr5_region is not None:
                ret.append(utr5_region)
        return ret

    def utr3_regions(self):
        """Return list of UTR3 regions

            Notes:
                Returned is a list of (x,y) pairs where coordinates x and y are 0-based, inclusive for x and exclusive for y
        """

        if self._any_unset(['exons', 'coding_end', 'strand']):
            return

        ret = []
        for exon in self.exons:
            utr3_region = exon.get_utr3(self.coding_end, self.strand)
            if utr3_region is not None:
                ret.append(utr3_region)
        return ret

    def get_cdna_length(self):
        """Return cDNA length"""

        if self.exons is None: return
        ret = 0
        for e in self.exons: ret += e.length()
        return ret

    def get_cds_length(self):
        """Return CDS length"""

        cds_regions = self.cds_regions()
        if cds_regions is None:
            return
        ret = 0
        for (cds_start, cds_end) in cds_regions:
            ret += cds_end - cds_start
        return ret

    def get_protein_length(self):
        """Return protein length"""

        cds_length = self.get_cds_length()
        if cds_length is None:
            return
        return int((cds_length - 3) / 3)

    def finalize(self):
        """Finalize transcript (set cDNA and protein length, start/end coordinates and info)"""

        # Calculate cDNA and protein length
        self.cdna_length = self.get_cdna_length()
        self.prot_length = self.get_protein_length()

        if self._any_unset(['cdna_length', 'prot_length', 'strand', 'exons']): return

        # Set start and end fields of transcript
        if self.strand == '+':
            self.start = self.exons[0].start
            self.end = self.exons[-1].end
        else:
            self.start = self.exons[-1].start
            self.end = self.exons[0].end

        # Set info field
        self.set_info()

    def _any_unset(self, fields):
        """Check if any of the required fields is None"""

        return any([getattr(self, f) is None for f in fields])

    def __str__(self):
        """String representation"""

        return self.id + '(' + self.chrom + ':' + str(self.start) + '-' + str(self.end) + ')[' + self.info + ']'

    __repr__ = __str__


class Exon(object):
    """Class for a single exon"""

    def __init__(self, s):
        """Constructor of Exon class

            Notes:
                s is a string of the syntax "start-end", where coordinates are 0-based, inclusive for start, exclusive for end
        """

        [s, e] = s.split('-')
        self.start = int(s)
        self.end = int(e)

    def length(self):
        """Return length of exon"""

        return self.end - self.start

    def get_cds(self, coding_start, coding_end):
        """Return CDS interval or None if there is no CDS in the exon"""

        coding_min = min(coding_start, coding_end)
        coding_max = max(coding_start, coding_end)

        if self.end - 1 < coding_min or self.start > coding_max:
            return

        cds_start = max(self.start, coding_min)
        cds_end = min(self.end - 1, coding_max) + 1
        return (cds_start, cds_end)

    def get_utr5(self, coding_start, strand):
        """Return UTR5 interval or None if there is no UTR5 in the exon"""

        if strand == '+':
            if self.start >= coding_start:
                return
            return (self.start, min(coding_start - 1, self.end - 1) + 1)
        else:
            if self.end - 1 <= coding_start:
                return
            return (max(self.start, coding_start + 1), self.end)

    def get_utr3(self, coding_end, strand):
        """Return UTR3 interval or None if there is no UTR3 in the exon"""

        if strand == '+':
            if self.end - 1 <= coding_end:
                return
            return (max(self.start, coding_end + 1), self.end)
        else:
            if self.start >= coding_end:
                return
            return (self.start, min(self.end - 1, coding_end - 1) + 1)

    def __str__(self):
        """String representation"""

        return 'exon:' + str(self.start) + '-' + str(self.end)

    __repr__ = __str__


class TranscriptDB(object):
    """Class for a transcript database"""

    # Allowed chromosome names
    allowed_chroms = map(str, range(1, 24)) + ['X', 'Y', 'MT']

    def __init__(self, filename):
        """Constructor of the TranscriptDB class"""

        self.source = ''
        self.date = ''
        self.build = ''
        self._filename = filename
        self._tabixfile = pysam.Tabixfile(filename)
        self._columns = []
        self._data = dict()
        self._read_header()

    def read(self):
        """Read transcript database from file"""

        for line in gzip.open(self._filename, 'r'):
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            cols = line.split('\t')
            self._data[cols[self._columns.index('ID')]] = line

    def contains(self, transcript_id):
        """Return True if transcript ID is found in the database"""

        return transcript_id in self._data

    def by_id(self, transcript_id):
        """Return transcript by ID"""

        ret = Transcript()
        ret.read_from_database_record(self._to_dict(self._data[transcript_id]))
        return ret

    def search_position(self, chrom, pos):
        """Search transcript database by genomic position"""

        ret = []
        reg = chrom + ':' + str(pos) + '-' + str(pos)
        lines = self._tabixfile.fetch(region=reg)
        for line in lines:
            t = Transcript()
            t.read_from_database_record(self._to_dict(line))
            ret.append(t)
        return ret

    def generator(self):
        """Return transcripts as a generator object"""

        for line in gzip.open(self._filename, 'r'):
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            ret = Transcript()
            ret.read_from_database_record(self._to_dict(line))
            yield ret

    def _read_header(self):
        """Read header information from file"""

        for line in gzip.open(self._filename, 'r'):
            line = line.strip()
            if line == '':
                continue
            if not line.startswith('#'):
                break
            if line.startswith('#source:'):
                self.source = line[9:].strip()
            elif line.startswith('#date:'):
                self.date = line[7:].strip()
            elif line.startswith('#build:'):
                self.build = line[8:].strip()
            else:
                self._columns = line[1:].split('\t')

    def _to_dict(self, line):
        """Convert line to dictionary"""

        ret = dict()
        cols = line.split('\t')
        for i in range(len(self._columns)):
            c = self._columns[i].lower()
            ret[c] = cols[i]
        return ret


class TranscriptDBWriter(object):
    """Class for creating new transcript database"""

    def __init__(self, fn, source='', build='', columns=[]):
        """Constructor of the TranscriptDBWriter class"""

        self._fn = fn
        self._source = source
        self._build = build
        self._columns = [x.lower() for x in columns]
        self._records = {c: [] for c in TranscriptDB.allowed_chroms}
        self.idx_chrom = self._columns.index('chrom')
        self.idx_start = self._columns.index('start')
        self.idx_end = self._columns.index('end')

    def add(self, transcript):
        """Add transcript to DB"""

        record = []
        for c in self._columns:
            if c in ['exons', 'cdna_exons']:
                record.append(','.join([str(e.start) + '-' + str(e.end) for e in getattr(transcript, c.lower())]))
            elif c in ['start', 'end', 'coding_start', 'coding_end', 'cdna_coding_start', 'cdna_coding_end']:
                record.append(int(getattr(transcript, c.lower())))
            else:
                record.append(str(getattr(transcript, c.lower())))
        self._records[transcript.chrom].append(record)

    def _sort_records(self):
        """Sort records by chrom, start, end"""

        idx_start = self._columns.index('start')
        idx_end = self._columns.index('end')
        for c in TranscriptDB.allowed_chroms:
            if c in self._records:
                self._records[c] = sorted(self._records[c], key=itemgetter(idx_start, idx_end))

    def _index_with_tabix(self):
        """Compress and index output file by Tabix"""

        pysam.tabix_compress(self._fn + '_tmp', self._fn + '.gz', force=True)
        pysam.tabix_index(self._fn + '.gz', seq_col=self.idx_chrom, start_col=self.idx_start, end_col=self.idx_end,
                          meta_char='#', force=True)

    def finalize(self):
        """Write to file, compress and index, clean up"""

        # Sort records by CHROM, START, END
        self._sort_records()

        # Initialize file and write header
        out = open(self._fn + '_tmp', 'w')
        out.write('#source: ' + self._source + '\n')
        out.write('#date: ' + str(datetime.datetime.today()).split()[0] + '\n')
        out.write('#build: ' + self._build + '\n')
        out.write('#' + '\t'.join([x.upper() for x in self._columns]) + '\n')

        # Write records to file
        for c in TranscriptDB.allowed_chroms:
            if c in self._records:
                for record in self._records[c]:
                    record = map(str, record)
                    out.write('\t'.join(record) + '\n')
        out.close()

        # Compress and index by Tabix
        self._index_with_tabix()

        # Remove temporary file
        os.remove(self._fn + '_tmp')

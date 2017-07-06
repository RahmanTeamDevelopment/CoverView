from collections import OrderedDict


class Transcript(object):
    def __init__(self, line):
        self.exons = []
        cols = line.split('\t')
        self.ENST = cols[0]
        self.geneSymbol = cols[1]
        self.geneID = cols[2]
        self.TRINFO = cols[3]
        self.chrom = cols[4]
        self.strand = int(cols[5])
        self.transcriptStart = int(cols[6])
        self.transcriptEnd = int(cols[7])
        self.codingStart = int(cols[8])
        self.codingStartGenomic = int(cols[9])
        self.codingEndGenomic = int(cols[10])
        # Initializing and adding exons
        for i in range(1, len(cols) - 11, 2):
            self.exons.append(Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i])))

    def getInfo(self):
        if self.strand == 1:
            ret = '+/'
        else:
            ret = '-/'
        cds = 0

        for exon in self.exons:
            cds += exon.length

        return ret + str(round((self.transcriptEnd - self.transcriptStart + 1) / 1000, 1)) + 'kb' + '/' + str(
            len(self.exons)) + '/' + str(round(cds / 1000, 1)) + 'kb'

    def isInUTR(self, pos):
        if self.strand == 1:
            return (pos < self.codingStartGenomic) or (pos > self.codingEndGenomic)
        else:
            return (pos > self.codingStartGenomic) or (pos < self.codingEndGenomic)


class Exon(object):
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start


def getTranscriptCoordinates(enstdb, chrom, pos):
    ret = OrderedDict()
    transcripts = findTranscripts(enstdb, chrom, pos)
    for enstid, transcript in transcripts.items():
        x, y = transformToCSNCoordinate(pos, transcript)
        transcoord = 'c.' + str(x)
        if y != 0:
            if y > 0:
                transcoord += '+' + str(y)
            else:
                transcoord += str(y)
        ret[transcript] = transcoord
    return ret


def findTranscripts(enstdb, chrom, pos):
    ret = OrderedDict()

    if chrom in enstdb.contigs:
        goodchrom = chrom
    else:
        if 'chr' + chrom in enstdb.contigs:
            goodchrom = 'chr' + chrom
        else:
            if chrom.startswith('chr') and chrom[3:] in enstdb.contigs:
                goodchrom = chrom[3:]
            else:
                return ret

    # Checking both end points of the variant
    reg = goodchrom + ':' + str(pos) + '-' + str(pos)
    hits = enstdb.fetch(region=reg)

    for line in hits:
        transcript = Transcript(line)
        if not (transcript.transcriptStart + 1 <= pos <= transcript.transcriptEnd):
            continue
        ret[transcript.ENST] = transcript

    return ret


# Transforming a genomic position ot CSN coordinate
def transformToCSNCoordinate(pos, transcript):
    prevExonEnd = 99999999
    # Checking if genomic position is within translated region
    if not transcript.isInUTR(pos):
        sumOfExonLengths = -transcript.codingStart + 1
        # Iterating through exons
        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]
            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prevExonEnd < pos < exon.start + 1:
                        if pos <= (exon.start + 1 - prevExonEnd) / 2 + prevExonEnd:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd
                        else:
                            x, y = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= (prevExonEnd - exon.end + 1) / 2 + exon.end:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos
                        else:
                            x, y = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos
            # Checking if genomic position is within exon
            if exon.start + 1 <= pos <= exon.end:
                if transcript.strand == 1:
                    return sumOfExonLengths + pos - exon.start, 0
                else:
                    return sumOfExonLengths + exon.end - pos + 1, 0
            # Calculating sum of exon lengths up to this point
            sumOfExonLengths += exon.length
            if transcript.strand == 1:
                prevExonEnd = exon.end
            else:
                prevExonEnd = exon.start + 1
    # If genomic position is within UTR
    else:
        if transcript.strand == 1:
            if pos < transcript.codingStartGenomic:
                return pos - transcript.codingStartGenomic, 0
            if pos > transcript.codingEndGenomic:
                return '+' + str(pos - transcript.codingEndGenomic), 0
        else:
            if pos > transcript.codingStartGenomic:
                return transcript.codingStartGenomic - pos, 0
            if pos < transcript.codingEndGenomic:
                return '+' + str(transcript.codingEndGenomic - pos), 0

import pysam

# Class representing the reference genome
class Reference(object):

    # Constructor
    def __init__(self, filename):
        # Openning tabix file representing the reference genome
        self.fastafile = pysam.Fastafile(filename)

    # Retrieving the sequence of a genomic region
    def getSequence(self, chrom, start, end):

        # Checking if chromosome name exists
        goodchrom = chrom
        if not goodchrom in self.fastafile.references:
            goodchrom = 'chr' + chrom
            if not goodchrom in self.fastafile.references:
                if chrom == 'MT':
                    goodchrom = 'chrM'
                    if not goodchrom in self.fastafile.references: return None
                else:
                    return None

        # Fetching data from reference genome
        if end < start: return ''
        if start < 1: start = 1

        if pysam.__version__ in ['0.7.7', '0.7.8', '0.8.0']:
            last = self.fastafile.getReferenceLength(goodchrom)
        else:
            last = self.fastafile.get_reference_length(goodchrom)

        if end > last: end = last
        seq = self.fastafile.fetch(goodchrom, start - 1, end)
        return seq.upper()


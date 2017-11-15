import interval


class BedFileParser(object):
    def __init__(self, bed_file):
        self.bed_file = bed_file

    def __iter__(self):
        return self

    def next(self):
        line = self.bed_file.next()
        cols = line.strip().split("\t")

        if len(cols) < 3:
            raise StandardError("Invalid line in BED file. Lines must have >= 3 columns")

        chromosome = cols[0]
        start_pos = int(cols[1])
        end_pos = int(cols[2])
        name = None

        if len(cols) > 3:
            name = cols[3]

        return interval.GenomicInterval(chromosome, start_pos, end_pos, name)

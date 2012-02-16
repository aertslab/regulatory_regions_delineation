#!/usr/bin/env python
import sys, re, operator


#col1:   seqid -> chrom
#col2:   if "protein_coding" continue
#       if "tRNA" skip line (only one line per tRNA)
#col3:   if "exon"       - count exons
#                       - record start/end (col4 & col5)
#                       - min start/end -> txStart; max start/end -> txEnd
#       if "CDS"        - min start/end -> cdsStart; max start/end -> cdsEnd
#col7:   strand
#col9:   attributes      format: key "value";
#                       gene_id -> ID


class Transcript(object):
    @staticmethod
    def _iterate_gtf(filename):
        with open(filename, 'r') as handle:
            for line in handle:
                if line.startswith('#'): continue
                seqid, source, type, start, end, score, strand, phase, attributes_str = line.rstrip().split('\t')
                if source != 'protein_coding': continue
                attributes = dict()
                for element in attributes_str.split(';'):
                    match = re.match("^\W*([A-Za-z_]+)\W+\"([A-Za-z_0-9\.]+)\"\W*$", element)
                    if match:
                        key, value = match.groups()
                        attributes[key] = value
                # Make intervals half-open and zero-based (in GFF they are closed intervals and 1-based) ...
                yield type, seqid, int(start)-1, int(end), strand, attributes

    @staticmethod
    def load_transcripts_from_file(filename):
        ID2Transcript = dict()
        for type, seqid, start, end, strand, attributes in Transcript._iterate_gtf(filename):
            ID = attributes['transcript_id']
            geneID = attributes['gene_id']
            transcript = ID2Transcript[ID] if ID in ID2Transcript else Transcript(ID, geneID, seqid, strand)
            if type == 'exon':
                exon = (start, end)
                transcript.exons.append(exon)
            elif type == 'CDS':
                CDS = (start, end)
                transcript.CDSs.append(CDS)
            if ID not in ID2Transcript: ID2Transcript[ID] = transcript
        return ID2Transcript.values()

    def __init__(self, ID, geneID, chromosome, strand):
        self.ID = ID
        self.geneID = geneID
        self.chromosome = chromosome
        self.strand = strand
        self.exons = []
        self.CDSs = []

    @property
    def txStart(self):
        if not self.exon_count: return None
        return min(self.exon_starts)

    @property
    def txEnd(self):
        if not self.exon_count: return None
        return max(self.exon_ends)

    @property
    def cdsStart(self):
        return min(map(operator.itemgetter(0), self.CDSs))

    @property
    def cdsEnd(self):
        return max(map(operator.itemgetter(1), self.CDSs))

    def __len__(self):
        if not self.exon_count: return 0
        return self.txEnd - self.txStart

    @property
    def exon_count(self):
        return len(self.exons)

    @property
    def exon_starts(self):
        if not self.exon_count: return tuple()
        return sorted(map(operator.itemgetter(0), self.exons))

    @property
    def exon_ends(self):
        if not self.exon_count: return tuple()
        return sorted(map(operator.itemgetter(1), self.exons))


def main():
    if len(sys.argv) != 2:
        print >>sys.stderr, "Wrong number of input arguments."
        sys.exit(2)
    filename = sys.argv[1]
    print >>sys.stderr, "Loading GFF into memory ..."
    transcripts = Transcript.load_transcripts_from_file(filename)
    print >>sys.stderr, "Conversion to tbl format ..."
    print >>sys.stdout, "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames"
    for tx in transcripts:
        if len(tx) == 0:
            print >>sys.stderr, "{0:s} is an empty transcript.".format(tx.ID)
            continue
        exonStarts = ','.join(map(str,tx.exon_starts)) + ','
        exonEnds = ','.join(map(str,tx.exon_ends)) + ','
        print >>sys.stdout, "NA\t{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t{5:d}\t{6:d}\t{7:d}\t{8:s}\t{9:s}\tNA\t{10:s}\tNA\tNA\tNA".format(
                tx.ID, tx.chromosome, tx.strand, tx.txStart, tx.txEnd, tx.cdsStart, tx.cdsEnd,
                tx.exon_count, exonStarts, exonEnds,
                tx.geneID)


if __name__ == "__main__":
    main()
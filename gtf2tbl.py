#!/usr/bin/env python
import operator
import re
import sys


# col1: seqid -> chrom
# col2: if "protein_coding" continue
#       if "tRNA" skip line (only one line per tRNA)
# col3: if "exon"       - count exons
#                       - record start/end (col4 & col5)
#                       - min start/end -> txStart; max start/end -> txEnd
#       if "CDS"        - min start/end -> cdsStart; max start/end -> cdsEnd
# col7: strand
# col9: attributes      format: key "value";
#                       gene_id -> ID


class Transcript(object):
    @staticmethod
    def _iterate_gtf(filename):
        with open(filename, 'r') as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                seq_id, source, type, start, end, score, strand, phase, attributes_str = line.rstrip().split('\t')
                if source != 'protein_coding':
                    continue
                attributes = dict()
                for element in attributes_str.split(';'):
                    match = re.match("^\W*([A-Za-z_]+)\W+\"([A-Za-z_0-9\.]+)\"\W*$", element)
                    if match:
                        key, value = match.groups()
                        attributes[key] = value
                # Make intervals half-open and zero-based (in GFF they are closed intervals and 1-based) ...
                yield type, seq_id, int(start) - 1, int(end), strand, attributes

    @staticmethod
    def load_transcripts_from_file(gtf_filename):
        transcript_id2Transcript = dict()
        for type, seq_id, start, end, strand, attributes in Transcript._iterate_gtf(gtf_filename):
            transcript_id = attributes['transcript_id']
            gene_id = attributes['gene_id']
            transcript = transcript_id2Transcript[
                transcript_id] if transcript_id in transcript_id2Transcript else Transcript(transcript_id, gene_id,
                                                                                            seq_id, strand)
            if type == 'exon':
                exon = (start, end)
                transcript.exons.append(exon)
            elif type == 'CDS':
                CDS = (start, end)
                transcript.CDSs.append(CDS)
            if transcript_id not in transcript_id2Transcript:
                transcript_id2Transcript[transcript_id] = transcript
        return transcript_id2Transcript.values()

    def __init__(self, transcript_id, gene_id, chromosome, strand):
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = []
        self.CDSs = []

    @property
    def tx_start(self):
        if not self.exon_count:
            return None
        return min(self.exon_starts)

    @property
    def tx_end(self):
        if not self.exon_count:
            return None
        return max(self.exon_ends)

    @property
    def cds_start(self):
        if not self.CDSs:
            return self.tx_start
        return min(map(operator.itemgetter(0), self.CDSs))

    @property
    def cds_end(self):
        if not self.CDSs:
            return self.tx_start
        return max(map(operator.itemgetter(1), self.CDSs))

    def __len__(self):
        if not self.exon_count:
            return 0
        return self.tx_end - self.tx_start

    @property
    def exon_count(self):
        return len(self.exons)

    @property
    def exon_starts(self):
        if not self.exon_count:
            return tuple()
        return sorted(map(operator.itemgetter(0), self.exons))

    @property
    def exon_ends(self):
        if not self.exon_count:
            return tuple()
        return sorted(map(operator.itemgetter(1), self.exons))


def main():
    if len(sys.argv) != 2:
        print('Wrong number of input arguments.', file=sys.stderr)
        sys.exit(2)
    gtf_filename = sys.argv[1]

    print('Loading GTF into memory ...', file=sys.stderr)
    transcripts = Transcript.load_transcripts_from_file(gtf_filename)

    print('Conversion to tbl format ...', file=sys.stderr)

    print('#bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
          'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames',
          sep='\t')

    for tx in transcripts:
        if len(tx) == 0:
            print('{0:s} is an empty transcript.'.format(tx.transcript_id), file=sys.stderr)
            continue

        exon_starts = ','.join(map(str, tx.exon_starts)) + ','
        exon_ends = ','.join(map(str, tx.exon_ends)) + ','

        print('NA', tx.transcript_id, tx.chromosome, tx.strand, tx.tx_start, tx.tx_end, tx.cds_start, tx.cds_end,
              tx.exon_count, exon_starts, exon_ends, 'NA', tx.gene_id, 'NA', 'NA', 'NA',
              sep='\t')


if __name__ == "__main__":
    main()

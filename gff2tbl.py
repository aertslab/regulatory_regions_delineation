#!/usr/bin/env python
import sys, operator


INCLUDE_TYPES = ('mRNA', 'miRNA', 'ncRNA', 'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'CDS', 'exon')


def iterate_gff_body(filename):
    with open(filename, 'r') as handle:
        for line in handle:
            if line.startswith('###'): return
            if line.startswith('##'): continue
            seqid, source, type, start, end, score, strand, phase, attributes_str = line.rstrip().split('\t')
            if type not in INCLUDE_TYPES: continue
            attributes = dict()
            for element in attributes_str.split(';'):
                key, value = element.split('=')
                if key == 'Parent': attributes[key] = value.split(',')
                elif key == 'ID': attributes[key] = value
                elif key == 'Name': attributes[key] = value
            # Make intervals half-open and zero-based (in GFF they are closed intervals and 1-based) ...
            yield type, seqid, int(start)-1, int(end), strand, attributes

    
class Entry:
    @staticmethod
    def from_mRNA(ID, chromosome, tx, strand, geneID):
        entry = Entry(ID)
        entry.chromosome, entry.tx, entry.strand, entry.geneID = chromosome, tx, strand, geneID
        return entry
    
    def __init__(self, ID):
        self.ID, self.chromosome, self.tx, self.strand, self.CDSs, self.exons, self.geneID = ID, '', (), '', [], [], ''
    
    def isempty(self):
        return len(self.tx) == 0
    
    def exon_count(self):
        return len(self.exons)
    
    def exon_starts(self):
        if self.exon_count() == 0: return tuple()
        return sorted(map(operator.itemgetter(0), self.exons))
        
    def exon_ends(self):
        if self.exon_count() == 0: return tuple()
        return sorted(map(operator.itemgetter(1), self.exons))
    
    def cds_span(self):
        if len(self.CDSs) == 0: return self.tx[1], self.tx[1]
        return min(map(operator.itemgetter(0), self.CDSs)), max(map(operator.itemgetter(1), self.CDSs))            


def load_data(filename):
    ID2entry = dict()
    for type, chromosome, start, end, strand, attributes in iterate_gff_body(filename):
        if type.endswith('RNA'):
            #Analyze attributes ...
            if 'ID' not in attributes:
                print >>sys.stderr, "'{0:s}' entry (ID: {1:s}) has no ID attribute.".format(type, ID)
                sys.exit(1)
            ID = attributes['ID']
            if 'Parent' in attributes and len(attributes['Parent']) == 1:
                geneID = attributes['Parent'][0]
            elif 'Name' in attributes: geneID = attributes['Name']
            else:    
                print >>sys.stderr, "'{0:s}' entry (ID: {1:s}) cannot be associated with a geneID.".format(type, ID)
                sys.exit(1)
            
            if ID in ID2entry.keys():
                #if ID2entry[ID].chromosome != chromosome or ID2entry[ID].strand != strand:
                #    print >>sys.stderr, "'mRNA' entry not compatible with child entry."
                #    sys.exit(1)
                ID2entry[ID].chromosome = chromosome
                ID2entry[ID].tx = (start, end)
                ID2entry[ID].strand = strand
                ID2entry[ID].geneID = geneID
            else: ID2entry[ID] = Entry.from_mRNA(ID, chromosome, (start, end), strand, geneID)
        elif type == 'exon':
            if 'Parent' not in attributes:
                print >>sys.stderr, "'exon' entry (ID: {0:s}) has no Parent attribute.".format(ID)
                sys.exit(1)
            exon = (start, end)
            for parentID in attributes['Parent']:
                if parentID not in ID2entry.keys(): ID2entry[parentID] = Entry(parentID)
                #elif ID2entry[parentID].chromosome != chromosome or ID2entry[parentID].strand != strand:
                #    print >>sys.stderr, "'exon' entry not compatible with parent 'mRNA' entry {0:s}.".format(parentID)
                #    sys.exit(1)
                ID2entry[parentID].exons.append(exon)
                #ID2entry[parentID].chromosome = chromosome
                #ID2entry[parentID].strand = strand
        elif type == 'CDS':
            if 'Parent' not in attributes:
                print >>sys.stderr, "'CDS' entry (ID: {0:s}) has no Parent attribute.".format(ID)
                sys.exit(1)
            CDS = (start, end)
            for parentID in attributes['Parent']:
                if parentID not in ID2entry.keys(): ID2entry[parentID] = Entry(parentID)
                #elif ID2entry[parentID].chromosome != chromosome or ID2entry[parentID].strand != strand:
                #    print >>sys.stderr, "'CDS' entry not compatible with parent 'mRNA' entry {0:s}.".format(parentID)
                #    sys.exit(1)
                ID2entry[parentID].CDSs.append(CDS)
                #ID2entry[parentID].chromosome = chromosome
                #ID2entry[parentID].strand = strand
    return ID2entry


def main():
    if len(sys.argv) != 2:
        print >>sys.stderr, "Wrong number of input arguments."
        sys.exit(2)
    filename = sys.argv[1]
    print >>sys.stderr, "Loading GFF into memory ..."
    ID2entry = load_data(filename)
    print >>sys.stderr, "Conversion to tbl format ..."
    print >>sys.stdout, "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames"
    for ID, entry in ID2entry.iteritems():
        if entry.isempty():
            print >>sys.stderr, "{0:s} is an empty entry.".format(ID)
            continue
        if entry.exon_count() == 0: entry.exons.append(entry.tx)
        cds = entry.cds_span()
        exonStarts = ','.join(map(str,entry.exon_starts())) + ','
        exonEnds = ','.join(map(str,entry.exon_ends())) + ','
        print >>sys.stdout, "NA\t{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t{5:d}\t{6:d}\t{7:d}\t{8:s}\t{9:s}\tNA\t{10:s}\tNA\tNA\tNA".format(
                ID, entry.chromosome, entry.strand, entry.tx[0], entry.tx[1], cds[0], cds[1],
                entry.exon_count(), exonStarts, exonEnds,
                entry.geneID)


if __name__ == "__main__":
    main()
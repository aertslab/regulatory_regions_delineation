#!/usr/bin/env python
import sys, os, sqlite3, operator
from transcript import Transcript
from interval import Interval
from piecewiselocation import PieceWiseLocation

# Attention! 08/06/2011: Different regulatory region delineation for
# non-coding sequences ...


LOAD_GENE_IDS_STATEMENT = r"""SELECT geneID FROM genes"""
GET_CHROMOSOME_LENGTHS_STATEMENT = r"""SELECT name, length FROM chromosomes"""

class Modes:
    FULL_TRANSCRIPT = 0
    ALL_INTRONS = 1
    UTR5_INTRON1 = 2
    NO_TRANSCRIPT = 3
    UTR5 = 4

        
def load_chromosome_lengths(connection):
    chromosome2length = dict()
    cursor = connection.cursor()
    cursor.execute(GET_CHROMOSOME_LENGTHS_STATEMENT)
    for name, length in cursor: chromosome2length[name] = length
    cursor.close()
    return chromosome2length


def chromosome_iterator(filename):
    with open(filename, 'r') as handle:
        for line in handle: yield line.rstrip()
    

def gene_ids_iterator(connection):
    cursor = connection.cursor()
    cursor.execute(LOAD_GENE_IDS_STATEMENT)
    for gene_id, in cursor: yield gene_id
    cursor.close()


def find_unique_transcript(transcripts, chromosomes):
    if len(transcripts) == 1: return transcripts[0]
    else:
        filtered = filter(lambda tx: tx.chromosome in chromosomes, transcripts)
        return filtered[0] if len(filtered) == 1 else None
    

def regulatory_regions_iterator(connection, chromosomes, chromosome2length,
            upstream_extension, downstream_extension, intronic_extension, mode):
    gene_ids = set(gene_ids_iterator(connection))
    gene_name2locations = dict()
    for gene_id in gene_ids:
        tx = find_unique_transcript(Transcript.load_by_gene_id(connection, gene_id), chromosomes)
        if not tx:
            print >>sys.stderr, "Skipped {0:s}: no unique transcript.".format(gene_id)
            continue
        
        if tx.coding_sequence().isempty():
            # Non-coding sequences = transcript is assumed not to be regulatory active
            # for the transcript itself...
            intragenic_location = tx.coding_sequence() # = empty!
        elif mode == Modes.FULL_TRANSCRIPT: # Intragenic location = full transcript ...
            intragenic_location = tx.transcript()
        elif mode == Modes.ALL_INTRONS: # Intragenic location = only introns ...
            introns = tx.introns_in_cds()
            if len(introns) > 0:
                intragenic_location = introns.location(0)
                for idx in range(1, len(introns)): intragenic_location += introns.location(idx)
            # Changed on 27/06/2011: 
            #else: intragenic_location = introns
            else: intragenic_location = tx.empty_interval()
        elif mode == Modes.UTR5_INTRON1: # Intragenic location = 5'UTR and 1st intron in CDS ...
            intragenic_location = tx.five_prime_utr()
            introns = tx.introns_in_cds()
            if len(introns) > 0:
                if tx.on_positive_strand(): intragenic_location += introns.location(0)
                else: intragenic_location += introns.location(len(introns) - 1)
        elif mode == Modes.UTR5: # Intragenic location = 5'UTR ...
            intragenic_location = tx.five_prime_utr()
        else: # No transcript ...
            intragenic_location = tx.empty_interval()

        if not intragenic_location.isempty() and intronic_extension > 0:
            if tx.on_positive_strand():
                filter_interval = Interval(tx.tss_as_bp_location(), tx.tss_as_bp_location() + intronic_extension)
                intragenic_location = intragenic_location.interval_limit(filter_interval)
            else:
                filter_interval = Interval(tx.tss_as_bp_location() - intronic_extension + 1, tx.tss_as_bp_location() + 1)
                intragenic_location = intragenic_location.interval_limit(filter_interval)
                
        # Remove coding sequences from all transcripts, including the current transcript
        # itself. The non-coding transcripts are not removed ...
        # Changed on 27/06/2011: 
        #if len(intragenic_location) > 0:
        if not intragenic_location.isempty():
            for rtx in Transcript.load_by_location(connection, intragenic_location.chromosome, intragenic_location.span()):
                intragenic_location -= rtx.coding_exons()
        regulatory_location = intragenic_location
        
        # Upstream region ...
        if upstream_extension > 0:
            upstream_location = tx.tss_shifted_1bp_upstream().extend_upstream(upstream_extension, chromosome2length)
            # Changed on 02/09/2011: check if upstream location is not empty ...
            if not upstream_location.isempty():
                for rtx in Transcript.load_by_location(connection, upstream_location.chromosome, upstream_location.span()):
                    upstream_location -= rtx.coding_exons() if tx.tss_as_bp_location() in rtx.transcript() else rtx.transcript()
                upstream_location = upstream_location.filter(tx.tss_shifted_1bp_upstream_as_bp_location())
                regulatory_location += upstream_location
        
        # Downstream region ...
        if downstream_extension > 0:
            downstream_location = tx.tes_shifted_1bp_downstream().extend_downstream(downstream_extension, chromosome2length)
            # Changed on 02/09/2011: check if downstream location is not empty ...
            if not downstream_location.isempty():
                for rtx in Transcript.load_by_location(connection, downstream_location.chromosome, downstream_location.span()):
                    downstream_location -= rtx.coding_exons() if tx.tes_as_bp_location() in rtx.transcript() else rtx.transcript()
                downstream_location = downstream_location.filter(tx.tes_shifted_1bp_downstream_as_bp_location())
                regulatory_location += downstream_location
        
        gene_name2locations.setdefault(tx.gene_name, []).append(regulatory_location)
    for gene_name in sorted(gene_name2locations.keys()):
        locations = gene_name2locations[gene_name]
        if (len(set([loc.chromosome for loc in locations])) != 1 or
            len(set([loc.on_positive_strand for loc in locations])) != 1):
            print >>sys.stderr, "Skipped {0:s}: cannot combine regulatory regions.".format(gene_name)
            continue
        combined_location = reduce(operator.add, locations)
        if combined_location.isempty():
            print >>sys.stderr, "Skipped {0:s}: no regulatory region remains.".format(gene_name)
        else:
            for idx, interval in enumerate(combined_location):
                yield (combined_location.chromosome, interval.start, interval.end,
                   "{0:s}#{1:d}".format(gene_name, idx + 1),
                   '+' if combined_location.on_positive_strand else '-')


def display_usage():
    print "Usage: python create-regulatory-regions-bed.py <sqlite3-db> <chromosomes> <upstream-extend> <downstream-extend> <intronic-extend> <full-transcript?>"


if __name__ == "__main__":
    if len(sys.argv) != 7:
        display_usage()
        print >>sys.stderr, "Wrong number of arguments."
        sys.exit(2)
        
    database_filename = sys.argv[1]
    if not os.path.isfile(database_filename):
        print >>sys.stderr, "Database file doesn't exist."
        sys.exit(1)
    chromosomes_filename = sys.argv[2]
    if not os.path.isfile(chromosomes_filename):
        print >>sys.stderr, "Chromosomes file doesn't exist."
        sys.exit(1)
    upstream_extend = int(sys.argv[3])
    downstream_extend = int(sys.argv[4])
    intronic_extend = int(sys.argv[5])
    if sys.argv[6] == "FullTx":
        mode = Modes.FULL_TRANSCRIPT
    elif sys.argv[6] == "AllIntrons":
        mode = Modes.ALL_INTRONS
    elif sys.argv[6] == "5utrIntron1":
        mode = Modes.UTR5_INTRON1
    elif sys.argv[6] == "NoTx":
        mode = Modes.NO_TRANSCRIPT
    elif sys.argv[6] == "5utr":
        mode = Modes.UTR5
    else:
        print >>sys.stderr, "'{0:s}' is an unknown mode.".format(sys.argv[6])
        sys.exit(1)
    
    chromosomes = set(chromosome_iterator(chromosomes_filename))
    with sqlite3.connect(database_filename) as connection:
        chromosome2length = load_chromosome_lengths(connection)
        for columns in regulatory_regions_iterator(connection, chromosomes, chromosome2length,
                    upstream_extend, downstream_extend, intronic_extend, mode):
            print >>sys.stdout, "{0:s}\t{1:d}\t{2:d}\t{3:s}\t0\t{4:s}".format(*columns)

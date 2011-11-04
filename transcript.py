import sqlite3
from interval import Interval
from piecewiselocation import PieceWiseLocation

EXON_SEPARATOR = ','
LOAD_TRANSCRIPTS_BY_LOCATION_STATEMENT = \
r"""SELECT * FROM genes WHERE chromosome = ? AND strand = ? AND txStart < ? AND txEnd > ?"""
LOAD_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT = \
r"""SELECT * FROM genes WHERE chromosome = ? AND txStart < ? AND txEnd > ?"""
LOAD_TRANSCRIPTS_BY_ID_STATEMENT = r"""SELECT * FROM genes WHERE ID = ?"""


class Transcript:
    # Changed on 27/06/2011: Creation of empty intervals ...
    def empty_interval(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart, self.txStart)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txEnd-1, self.txEnd-1)
    
    @staticmethod
    def load_by_id(connection, id):
        cursor = connection.cursor()
        cursor.execute(LOAD_TRANSCRIPTS_BY_ID_STATEMENT, (id,))
        transcripts = [Transcript(*columns) for columns in cursor]
        cursor.close()
        return transcripts


    @staticmethod
    def load_by_location(connection, chromosome, interval, strand=None):
        cursor = connection.cursor()
        if strand:
            #CAVE: end and start should be supplied to the statement in opposite order ...
            cursor.execute(LOAD_TRANSCRIPTS_BY_LOCATION_STATEMENT, (chromosome, strand, interval.end, interval.start)) 
        else:
            #CAVE: end and start should be supplied to the statement in opposite order ...
            cursor.execute(LOAD_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT, (chromosome, interval.end, interval.start))
        transcripts = [Transcript(*columns) for columns in cursor]
        cursor.close()
        return transcripts

    
    def __init__(self, ID, chromosome, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, altName):
        self.ID = ID
        self.altName = altName
        self.chromosome = chromosome
        self.strand = strand
        self.txStart = txStart
        self.txEnd = txEnd
        self.cdsStart = cdsStart
        self.cdsEnd = cdsEnd
        self.exonCount = exonCount
        self._exonStarts = exonStarts
        self._exonEnds = exonEnds

    
    def __str__(self):
        return "{0:s}\t{1:d}\t{2:d}\t{3:s}\t{4:s}".format(self.chromosome, self.txStart, self.txEnd, self.ID, self.strand)


    def on_positive_strand(self):
        return self.strand == '+'
        
        
    def on_negative_strand(self):
        return self.strand == '-'
        
    
    def tss(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart, self.txStart+1)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txEnd - 1, self.txEnd)
    
    
    def tss_as_bp_location(self):
        return self.txStart if self.on_positive_strand() else self.txEnd - 1
    
    
    def tes(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txEnd - 1, self.txEnd)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart, self.txStart + 1)
    
    
    def tes_as_bp_location(self):
        return self.txEnd - 1 if self.on_positive_strand() else self.txStart
    
    
    def tss_shifted_1bp_upstream(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart - 1, self.txStart)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txEnd, self.txEnd + 1)
    
    
    def tss_shifted_1bp_upstream_as_bp_location(self):
        return self.txStart - 1 if self.on_positive_strand() else self.txEnd
    
    
    def tes_shifted_1bp_downstream(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txEnd, self.txEnd + 1)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart - 1, self.txStart)
    
    
    def tes_shifted_1bp_downstream_as_bp_location(self):
        return self.txEnd if self.on_positive_strand() else self.txStart - 1
    
    
    def transcript(self):
        return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart, self.txEnd)
    
    
    def five_prime_utr(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart, self.cdsStart)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.cdsEnd, self.txEnd)


    def three_prime_utr(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.cdsEnd, self.txEnd)
        else:
            return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.txStart, self.cdsStart)

    
    def coding_sequence(self):
        return PieceWiseLocation.singleton(self.chromosome, self.on_positive_strand(), self.cdsStart, self.cdsEnd)

    
    def _introns_interval_iterator(self):
        #CAVE: list of exon start and ends has an additional final comma ...
        for interval in zip(map(int, self._exonEnds.split(EXON_SEPARATOR)[0:-2]),
                            map(int, self._exonStarts.split(EXON_SEPARATOR)[1:-1])):
            yield interval


    def introns(self):
        introns = [Interval(start, end) for start, end in self._introns_interval_iterator()]
        return PieceWiseLocation(self.chromosome, self.on_positive_strand(), introns)
        

    def introns_in_cds(self):
        introns = [Interval(max(intronStart, self.cdsStart), min(intronEnd, self.cdsEnd))
                 for intronStart, intronEnd in self._introns_interval_iterator()
                 if intronEnd > self.cdsStart and intronStart < self.cdsEnd]
        return PieceWiseLocation(self.chromosome, self.on_positive_strand(), introns)
    
    
    def _exons_interval_iterator(self):
        #CAVE: list of exon start and ends has an additional final comma ...
        for interval in zip(map(int, self._exonStarts.split(EXON_SEPARATOR)[0:-1]),
                              map(int, self._exonEnds.split(EXON_SEPARATOR)[0:-1])):
            yield interval
    
    
    def exons(self):
        exons = [Interval(start, end) for start, end in self._exons_interval_iterator()]
        return PieceWiseLocation(self.chromosome, self.on_positive_strand(), exons)
        
    
    def coding_exons(self):
        exons = [Interval(max(exonStart, self.cdsStart), min(exonEnd, self.cdsEnd))
                 for exonStart, exonEnd in self._exons_interval_iterator()
                 if exonEnd > self.cdsStart and exonStart < self.cdsEnd]
        return PieceWiseLocation(self.chromosome, self.on_positive_strand(), exons)
from interval import Interval
from piecewiselocation import PieceWiseLocation

EXON_SEPARATOR = ","

LOAD_ALL_TRANSCRIPTS_BY_LOCATION_STATEMENT = r"""
    SELECT * FROM genes
    WHERE chromosome = ? AND strand = ? AND txStart < ? AND txEnd > ?;
"""

LOAD_ONLY_NM_TRANSCRIPTS_BY_LOCATION_STATEMENT = r"""
    SELECT * FROM genes
    WHERE chromosome = ? AND strand = ? AND txStart < ? AND txEnd > ? AND geneID LIKE "NM_%";
"""

LOAD_ALL_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT = r"""
    SELECT * FROM genes
    WHERE chromosome = ? AND txStart < ? AND txEnd > ?;
"""

LOAD_ONLY_NM_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT = r"""
    SELECT * FROM genes
    WHERE chromosome = ? AND txStart < ? AND txEnd > ? AND geneID LIKE "NM_%";
"""

LOAD_TRANSCRIPTS_BY_GENE_ID_STATEMENT = r"""SELECT * FROM genes WHERE geneID = ?;"""


class Transcript:
    # Changed on 27/06/2011: Creation of empty intervals ...
    def empty_interval(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start,
                self.tx_start,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_end - 1,
                self.tx_end - 1,
            )

    @staticmethod
    def load_by_gene_id(connection, gene_id):
        cursor = connection.cursor()
        cursor.execute(LOAD_TRANSCRIPTS_BY_GENE_ID_STATEMENT, (gene_id,))
        transcripts = [Transcript(*columns) for columns in cursor]
        cursor.close()
        return transcripts

    @staticmethod
    def load_by_location(connection, chromosome, interval, strand=None, only_NMs=False):
        cursor = connection.cursor()

        if strand:
            LOAD_TRANSCRIPTS_BY_LOCATION_STATEMENT = (
                LOAD_ONLY_NM_TRANSCRIPTS_BY_LOCATION_STATEMENT
                if only_NMs
                else LOAD_ALL_TRANSCRIPTS_BY_LOCATION_STATEMENT
            )

            # CAVE: end and start should be supplied to the statement in opposite order ...
            cursor.execute(
                LOAD_TRANSCRIPTS_BY_LOCATION_STATEMENT,
                (chromosome, strand, interval.end, interval.start),
            )
        else:
            LOAD_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT = (
                LOAD_ONLY_NM_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT
                if only_NMs
                else LOAD_ALL_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT
            )

            # CAVE: end and start should be supplied to the statement in opposite order ...
            cursor.execute(
                LOAD_TRANSCRIPTS_BY_STRANDLESS_LOCATION_STATEMENT,
                (chromosome, interval.end, interval.start),
            )

        transcripts = [Transcript(*columns) for columns in cursor]
        cursor.close()
        return transcripts

    def __init__(
        self,
        gene_id,
        chromosome,
        strand,
        tx_start,
        tx_end,
        cds_start,
        cds_end,
        exon_count,
        exon_starts,
        exon_ends,
        gene_name,
    ):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chromosome = chromosome
        self.strand = strand
        self.tx_start = tx_start
        self.tx_end = tx_end
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.exon_count = exon_count
        self._exon_starts = exon_starts
        self._exon_ends = exon_ends

    def __str__(self):
        return "{0:s}\t{1:d}\t{2:d}\t{3:s}\t{4:s}".format(
            self.chromosome,
            self.tx_start,
            self.tx_end,
            self.gene_id,
            self.strand,
        )

    def on_positive_strand(self):
        return self.strand == "+"

    def on_negative_strand(self):
        return self.strand == "-"

    def cds_start_location(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.cds_start,
                self.cds_start + 1,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.cds_end - 1,
                self.cds_end,
            )

    def cds_start_as_bp_location(self):
        return self.cds_start if self.on_positive_strand() else self.cds_end - 1

    def tss(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start,
                self.tx_start + 1,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_end - 1,
                self.tx_end,
            )

    def tss_as_bp_location(self):
        return self.tx_start if self.on_positive_strand() else self.tx_end - 1

    def tes(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_end - 1,
                self.tx_end,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start,
                self.tx_start + 1,
            )

    def tes_as_bp_location(self):
        return self.tx_end - 1 if self.on_positive_strand() else self.tx_start

    def cds_start_shifted_1bp_upstream(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.cds_start - 1,
                self.cds_start,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.cds_end,
                self.cds_end + 1,
            )

    def cds_start_shifted_1bp_upstream_as_bp_location(self):
        return self.cds_start - 1 if self.on_positive_strand() else self.cds_end

    def tss_shifted_1bp_upstream(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start - 1,
                self.tx_start,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_end,
                self.tx_end + 1,
            )

    def tss_shifted_1bp_upstream_as_bp_location(self):
        return self.tx_start - 1 if self.on_positive_strand() else self.tx_end

    def tes_shifted_1bp_downstream(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_end,
                self.tx_end + 1,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start - 1,
                self.tx_start,
            )

    def tes_shifted_1bp_downstream_as_bp_location(self):
        return self.tx_end if self.on_positive_strand() else self.tx_start - 1

    def transcript(self):
        return PieceWiseLocation.singleton(
            self.chromosome,
            self.on_positive_strand(),
            self.tx_start,
            self.tx_end,
        )

    def five_prime_utr(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start,
                self.cds_start,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.cds_end,
                self.tx_end,
            )

    def three_prime_utr(self):
        if self.on_positive_strand():
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.cds_end,
                self.tx_end,
            )
        else:
            return PieceWiseLocation.singleton(
                self.chromosome,
                self.on_positive_strand(),
                self.tx_start,
                self.cds_start,
            )

    def coding_sequence(self):
        return PieceWiseLocation.singleton(
            self.chromosome,
            self.on_positive_strand(),
            self.cds_start,
            self.cds_end,
        )

    def _introns_interval_iterator(self):
        # CAVE: list of exon start and ends has an additional final comma ...
        for interval in zip(
            map(int, self._exon_ends.split(EXON_SEPARATOR)[0:-2]),
            map(int, self._exon_starts.split(EXON_SEPARATOR)[1:-1]),
        ):
            yield interval

    def introns(self):
        introns = [
            Interval(start, end) for start, end in self._introns_interval_iterator()
        ]
        return PieceWiseLocation(
            self.chromosome,
            self.on_positive_strand(),
            introns,
        )

    def introns_in_cds(self):
        introns = [
            Interval(max(intronStart, self.cds_start), min(intronEnd, self.cds_end))
            for intronStart, intronEnd in self._introns_interval_iterator()
            if intronEnd > self.cds_start and intronStart < self.cds_end
        ]
        return PieceWiseLocation(
            self.chromosome,
            self.on_positive_strand(),
            introns,
        )

    def _exons_interval_iterator(self):
        # CAVE: list of exon start and ends has an additional final comma ...
        for interval in zip(
            map(int, self._exon_starts.split(EXON_SEPARATOR)[0:-1]),
            map(int, self._exon_ends.split(EXON_SEPARATOR)[0:-1]),
        ):
            yield interval

    def exons(self):
        exons = [Interval(start, end) for start, end in self._exons_interval_iterator()]
        return PieceWiseLocation(
            self.chromosome,
            self.on_positive_strand(),
            exons,
        )

    def coding_exons(self):
        exons = [
            Interval(max(exonStart, self.cds_start), min(exonEnd, self.cds_end))
            for exonStart, exonEnd in self._exons_interval_iterator()
            if exonEnd > self.cds_start and exonStart < self.cds_end
        ]
        return PieceWiseLocation(
            self.chromosome,
            self.on_positive_strand(),
            exons,
        )

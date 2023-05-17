import sqlite3
import unittest

from interval import Interval
from piecewiselocation import PieceWiseLocation
from transcript import Transcript

DB_FILENAME = "hg19-refseq-genes.sqlite3.db"


class IntervalTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_contains(self):
        assert 1 in Interval(1, 3)
        assert 2 in Interval(1, 3)
        assert 3 not in Interval(1, 3)

    def test_eq(self):
        assert Interval(1, 3) == Interval(1, 3)
        assert Interval(1, 2) != Interval(1, 3)
        assert Interval(2, 3) != Interval(1, 3)
        assert Interval(1, 3) != Interval(4, 5)

    def test_isempty(self):
        assert Interval(1, 1).isempty()
        assert Interval(2, 1).isempty()
        assert not Interval(1, 2).isempty()
        assert not Interval(1, 5).isempty()

    def test_overlaps_with(self):
        assert Interval(1, 3).overlaps_with(Interval(2, 4))
        assert Interval(1, 5).overlaps_with(Interval(2, 4))
        assert Interval(1, 5).overlaps_with(Interval(1, 2))
        assert Interval(1, 5).overlaps_with(Interval(4, 5))
        assert not Interval(1, 5).overlaps_with(Interval(1, 1))
        assert not Interval(1, 5).overlaps_with(Interval(7, 7))
        assert not Interval(1, 3).overlaps_with(Interval(3, 5))

    def test_cmp(self):
        assert Interval(1, 2) < Interval(2, 5)
        assert Interval(1, 5) < Interval(1, 7)
        assert Interval(1, 5) <= Interval(1, 5)
        assert not Interval(1, 5) < Interval(1, 5)

    def test_merge(self):
        intervals = list(
            reversed(
                [
                    Interval(1, 3),
                    Interval(23, 23),
                    Interval(2, 5),
                    Interval(4, 7),
                    Interval(9, 10),
                    Interval(10, 11),
                    Interval(13, 16),
                    Interval(15, 17),
                ]
            )
        )
        result = Interval.merge(intervals)
        assert len(result) == 3
        assert result[0] == Interval(1, 7)
        assert result[1] == Interval(9, 11)
        assert result[2] == Interval(13, 17)
        result = Interval.merge([Interval(1, 2)])
        assert len(result) == 1
        assert result[0] == Interval(1, 2)


class PieceWiseLocationTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_span(self):
        assert Interval(1, 67) == PieceWiseLocation("chr1", True, [Interval(23, 67), Interval(14, 24), Interval(1, 4)]).span()

    def test_add(self):
        result = PieceWiseLocation(
            "chr1", True, [Interval(1, 5), Interval(10, 24), Interval(78, 100)]
        ) + PieceWiseLocation("chr1", True, [Interval(5, 6), Interval(26, 30)])
        assert len(result) == 4
        assert result[0] == Interval(1, 6)
        assert result[1] == Interval(10, 24)
        assert result[2] == Interval(26, 30)
        assert result[3] == Interval(78, 100)

    def test_sub(self):
        result = PieceWiseLocation(
            "chr1", True, [Interval(1, 4), Interval(10, 24), Interval(78, 100)]
        ) - PieceWiseLocation("chr1", True, [Interval(4, 5), Interval(26, 30)])
        assert len(result) == 3
        assert result[0] == Interval(1, 4)
        assert result[1] == Interval(10, 24)
        assert result[2] == Interval(78, 100)
        result = PieceWiseLocation(
            "chr1", True, [Interval(1, 8), Interval(10, 24), Interval(78, 100)]
        ) - PieceWiseLocation("chr1", True, [Interval(4, 5), Interval(26, 30)])
        assert len(result) == 4
        assert result[0] == Interval(1, 4)
        assert result[1] == Interval(5, 8)
        assert result[2] == Interval(10, 24)
        assert result[3] == Interval(78, 100)
        result = PieceWiseLocation(
            "chr1", True, [Interval(1, 6), Interval(10, 24), Interval(78, 100)]
        ) - PieceWiseLocation("chr1", True, [Interval(4, 13), Interval(21, 79)])
        assert len(result) == 3
        assert result[0] == Interval(1, 4)
        assert result[1] == Interval(13, 21)
        assert result[2] == Interval(79, 100)

    def test_extend_upstream(self):
        chromosome2length = {"chr1": 500}
        result = PieceWiseLocation(
            "chr1", True, [Interval(100, 200), Interval(300, 400)]
        ).extend_upstream(75, chromosome2length)
        assert len(result) == 2
        assert result[0] == Interval(25, 200)
        assert result[1] == Interval(300, 400)
        result = PieceWiseLocation(
            "chr1", True, [Interval(100, 200), Interval(300, 400)]
        ).extend_upstream(175, chromosome2length)
        assert len(result) == 2
        assert result[0] == Interval(0, 200)
        assert result[1] == Interval(300, 400)
        result = PieceWiseLocation(
            "chr1", False, [Interval(100, 200), Interval(300, 400)]
        ).extend_upstream(75, chromosome2length)
        assert len(result) == 2
        assert result[0] == Interval(100, 200)
        assert result[1] == Interval(300, 475)
        result = PieceWiseLocation(
            "chr1", False, [Interval(100, 200), Interval(300, 400)]
        ).extend_upstream(175, chromosome2length)
        assert len(result) == 2
        assert result[0] == Interval(100, 200)
        assert result[1] == Interval(300, 500)


class TranscriptTest(unittest.TestCase):
    def setUp(self):
        self.connection = sqlite3.connect(DB_FILENAME)

    def tearDown(self):
        self.connection.close()

    def test_NM_005649(self):
        # NM_005649    chr5    -    178138521    178157703    178139060    178156023
        # 5    178138521,178152376,178153999,178155990,178157556,    178140622,178152472,178154126,178156074,178157703,    ZNF354A
        tx = Transcript.load_by_gene_id(self.connection, "NM_005649")[0]
        assert tx.transcript()[0] == Interval(178138521, 178157703)
        assert tx.coding_sequence()[0] == Interval(178139060, 178156023)
        assert tx.five_prime_utr()[0] == Interval(178156023, 178157703)
        assert tx.three_prime_utr()[0] == Interval(178138521, 178139060)
        exons = tx.exons()
        assert len(exons) == 5
        assert exons[0] == Interval(178138521, 178140622)
        assert exons[4] == Interval(178157556, 178157703)
        introns = tx.introns()
        assert len(introns) == 4
        assert introns[0] == Interval(178140622, 178152376)
        assert introns[3] == Interval(178156074, 178157556)


if __name__ == "__main__":
    unittest.main()

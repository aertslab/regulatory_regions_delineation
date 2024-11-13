#!/usr/bin/env python3
import argparse
import operator
import re
import sys


# Convert GTF file to GenePred format.
#
# Rules for determining coordinates should be following the behaviour of
# Kent Tools gtfToGenePred:
#
#     gtfToGenePred -genePredExt genes.gtf genes.gp
#
#   - Use BED-like coordinates (zero-based half open intervals).
#   - Take protein coding transcripts only.
#   - If transcript does not have a CDS, set highest genomic coordinate
#     of exon as start and end for cdsStart and cdsEnd.
#   - Stop codon is included in cdsStart and cdsEnd. (CDS normally does
#     not include stop codon.)


class Transcript:
    @staticmethod
    def _iterate_gtf(filename):
        with open(filename, "r") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                (
                    seq_id,
                    source,
                    type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attributes_str,
                ) = line.rstrip().split("\t")
                attributes = {}
                for element in attributes_str.split(";"):
                    match = re.match(r'^\W*([A-Za-z_]+)\W+"([^"]+)"\W*$', element)
                    if match:
                        key, value = match.groups()
                        attributes[key] = value
                gene_biotype = attributes.get("gene_biotype")
                if gene_biotype and gene_biotype == "protein_coding":
                    # Make intervals half-open and zero-based (in GTF they are closed intervals and 1-based) ...
                    yield type, seq_id, int(start) - 1, int(end), strand, attributes

    @staticmethod
    def load_transcripts_from_file(gtf_filename):
        transcript_id2Transcript = {}
        for type, seq_id, start, end, strand, attributes in Transcript._iterate_gtf(
            gtf_filename
        ):
            if (
                type == "transcript"
                or type == "exon"
                or type == "CDS"
                or type == "stop_codon"
            ):
                transcript_id = attributes["transcript_id"]
                gene_id = attributes["gene_id"]
                transcript = (
                    transcript_id2Transcript[transcript_id]
                    if transcript_id in transcript_id2Transcript
                    else Transcript(transcript_id, gene_id, seq_id, strand)
                )
                if type == "exon":
                    exon = (start, end)
                    transcript.exons.append(exon)
                elif type == "CDS":
                    CDS = (start, end)
                    transcript.CDSs.append(CDS)
                elif type == "stop_codon":
                    stop_codon = (start, end)
                    transcript.stop_codon = stop_codon
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
        self.stop_codon = None

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
            # Kent Tools gtfToGenePred seems to put CDS start always a the highest
            # genomic coordinate of the transcript.
            return self.tx_end  # if self.strand == "-" else self.tx_start
        cds_start = min(map(operator.itemgetter(0), self.CDSs))

        if self.stop_codon and self.strand == "-":
            # If transcript has an annotated stop codon, set start position as CDS
            # start if transcript is on negative strand.
            return self.stop_codon[0]

        return cds_start

    @property
    def cds_end(self):
        if not self.CDSs:
            # Kent Tools gtfToGenePred seems to put CDS end always a the highest
            # genomic coordinate of the transcript.
            return self.tx_end  # if self.strand == "-" else self.tx_start
        cds_end = max(map(operator.itemgetter(1), self.CDSs))

        if self.stop_codon and self.strand == "+":
            # If transcript has an annotated stop codon, set end position as CDS
            # end if transcript is on positive strand.
            return self.stop_codon[1]
        return cds_end

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
            return ()
        return sorted(map(operator.itemgetter(0), self.exons))

    @property
    def exon_ends(self):
        if not self.exon_count:
            return ()
        return sorted(map(operator.itemgetter(1), self.exons))

    @property
    def exon_frames(self):
        """
        Get exon frames.

        Get for each exon how many nucleotides it should get from the previous
        exon to be in the open reading frame. If an exon is not in the CDS, "-1"
        is used.
        """
        if not self.exon_count:
            return ()
        exon_frames_list = []
        left_over_codon_nucleotides = 0
        if self.strand == "-":
            # Negative strand.
            for exon_start, exon_end in zip(
                self.exon_starts[::-1], self.exon_ends[::-1]
            ):
                if (
                    not self.CDSs
                    or (exon_start > self.cds_end)
                    or (exon_end < self.cds_start)
                ):
                    # Set exon frame for current exon to -1 if transcript does not
                    # have a CDS or exon is located outside of CDS.
                    exon_frames_list.append(-1)
                else:
                    exon_frames_list.append(left_over_codon_nucleotides)
                    if exon_end > self.cds_end:
                        exon_end = self.cds_end
                    # Calculate number of nucleotides that are left over in the
                    # current exon and that will be used in the next exon to continue
                    # the open reading frame.
                    left_over_codon_nucleotides = (
                        exon_end + left_over_codon_nucleotides - exon_start
                    ) % 3
            return exon_frames_list[::-1]
        else:
            # Positive strand.
            for exon_start, exon_end in zip(self.exon_starts, self.exon_ends):
                if (
                    not self.CDSs
                    or (exon_end < self.cds_start)
                    or (exon_start > self.cds_end)
                ):
                    # Set exon frame for current exon to -1 if transcript does not
                    # have a CDS or exon is located outside of CDS.
                    exon_frames_list.append(-1)
                else:
                    exon_frames_list.append(left_over_codon_nucleotides)
                    if exon_start < self.cds_start:
                        exon_start = self.cds_start
                    # Calculate number of nucleotides that are left over in the
                    # current exon and that will be used in the next exon to continue
                    # the open reading frame.
                    left_over_codon_nucleotides = (
                        exon_end + left_over_codon_nucleotides - exon_start
                    ) % 3
            return exon_frames_list


def main():
    parser = argparse.ArgumentParser(
        description="Convert a GTF file to a GenePred file."
    )

    parser.add_argument(
        "-i",
        "--gtf",
        dest="gtf_filename",
        action="store",
        type=str,
        required=True,
        help="GTF input filename.",
    )
    parser.add_argument(
        "-o",
        "--gp",
        dest="genepred_filename",
        action="store",
        type=str,
        required=True,
        help="GenePred (gene prediction) output filename.",
    )

    args = parser.parse_args()

    print("Loading GTF into memory ...", file=sys.stderr)
    transcripts = Transcript.load_transcripts_from_file(args.gtf_filename)

    print("Conversion to GenePred format ...", file=sys.stderr)

    with open(args.genepred_filename, "w") as fh_genepred:
        print(
            "#bin",
            "name",
            "chrom",
            "strand",
            "txStart",
            "txEnd",
            "cdsStart",
            "cdsEnd",
            "exonCount",
            "exonStarts",
            "exonEnds",
            "score",
            "name2",
            "cdsStartStat",
            "cdsEndStat",
            "exonFrames",
            sep="\t",
            file=fh_genepred,
        )

        for tx in transcripts:
            if len(tx) == 0:
                print(
                    f"{tx.transcript_id:s} is an empty transcript.",
                    file=sys.stderr,
                )
                continue

            exon_starts = ",".join(map(str, tx.exon_starts)) + ","
            exon_ends = ",".join(map(str, tx.exon_ends)) + ","
            exon_frames = ",".join(map(str, tx.exon_frames)) + ","

            print(
                "NA",
                tx.transcript_id,
                tx.chromosome,
                tx.strand,
                tx.tx_start,
                tx.tx_end,
                tx.cds_start,
                tx.cds_end,
                tx.exon_count,
                exon_starts,
                exon_ends,
                "NA",
                tx.gene_id,
                "NA",
                "NA",
                exon_frames,
                sep="\t",
                file=fh_genepred,
            )


if __name__ == "__main__":
    main()

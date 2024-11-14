#!/usr/bin/env python3
import argparse
import os
import sqlite3
import sys

from interval import Interval
from transcript import Transcript

# Attention! 08/06/2011: Different regulatory region delineation for non-coding sequences ...


LOAD_GENE_IDS_STATEMENT = r"""SELECT geneID FROM genes"""
GET_CHROMOSOME_LENGTHS_STATEMENT = r"""SELECT name, length FROM chromosomes"""


class Modes:
    FULL_TRANSCRIPT = 0
    ALL_INTRONS = 1
    UTR5_INTRON1 = 2
    NO_TRANSCRIPT = 3
    UTR5 = 4


def load_chromosome_lengths(connection):
    chromosome2length = {}
    cursor = connection.cursor()
    cursor.execute(GET_CHROMOSOME_LENGTHS_STATEMENT)
    for name, length in cursor:
        chromosome2length[name] = length
    cursor.close()
    return chromosome2length


def chromosome_iterator(filename, exclude_alt=True, exclude_random=True):
    with open(filename, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                # Skip comment lines.
                continue

            line = line.rstrip()

            if line:
                # Get first column as chromosome name.
                chrom = line.split("\t", 1)[0]

                # Skip alternative chromosomes, if requested.
                if exclude_alt and chrom.endswith("_alt"):
                    continue

                # Skip random chromosomes, if requested.
                if exclude_random and chrom.endswith("random"):
                    continue

                yield chrom


def gene_ids_iterator(connection):
    cursor = connection.cursor()
    cursor.execute(LOAD_GENE_IDS_STATEMENT)
    for (gene_id,) in cursor:
        yield gene_id
    cursor.close()


def find_unique_transcript(transcripts, chromosomes):
    transcripts = list(filter(lambda tx: tx.chromosome in chromosomes, transcripts))
    return transcripts[0] if len(transcripts) == 1 else None


def regulatory_regions_iterator(
    connection,
    chromosomes,
    chromosome2length,
    upstream_extension,
    downstream_extension,
    intronic_extension,
    mode,
):
    gene_ids = set(gene_ids_iterator(connection))
    gene_name2locations = {}
    for gene_id in gene_ids:
        # If gene_id is "NM_*", do not restrict regulatory locations by nearby "NR_*" but only by other "NM_*".
        is_NM = gene_id.startswith("NM_")

        tx = find_unique_transcript(
            Transcript.load_by_gene_id(connection, gene_id), chromosomes
        )
        if not tx:
            print(f"Skipped {gene_id:s}: no unique transcript.", file=sys.stderr)
            continue

        if mode == Modes.FULL_TRANSCRIPT:
            # Intragenic location = full transcript ...
            intragenic_location = tx.transcript()
        elif mode == Modes.ALL_INTRONS:
            # Intragenic location = only introns ...
            introns = tx.introns_in_cds()
            if len(introns) > 0:
                intragenic_location = introns.location(0)
                for idx in range(1, len(introns)):
                    intragenic_location += introns.location(idx)
            else:
                intragenic_location = tx.empty_interval()
        elif mode == Modes.UTR5_INTRON1:
            # Intragenic location = 5'UTR and 1st intron in CDS ...
            intragenic_location = tx.five_prime_utr()
            introns = tx.introns_in_cds()
            if len(introns) > 0:
                if tx.on_positive_strand():
                    intragenic_location += introns.location(0)
                else:
                    intragenic_location += introns.location(len(introns) - 1)
        elif mode == Modes.UTR5:
            # Intragenic location = 5'UTR ...
            intragenic_location = tx.five_prime_utr()
        else:
            # No transcript ...
            intragenic_location = tx.empty_interval()

        if not intragenic_location.isempty() and intronic_extension > 0:
            if tx.on_positive_strand():
                filter_interval = Interval(
                    tx.tss_as_bp_location(),
                    tx.tss_as_bp_location() + intronic_extension,
                )
                intragenic_location = intragenic_location.interval_limit(
                    filter_interval
                )
            else:
                filter_interval = Interval(
                    tx.tss_as_bp_location() - intronic_extension + 1,
                    tx.tss_as_bp_location() + 1,
                )
                intragenic_location = intragenic_location.interval_limit(
                    filter_interval
                )

        # Remove coding sequences from all transcripts, including the current transcript
        # itself. The non-coding transcripts are not removed ...
        if not intragenic_location.isempty():
            for rtx in Transcript.load_by_location(
                connection,
                intragenic_location.chromosome,
                intragenic_location.span(),
                only_NMs=is_NM,
            ):
                intragenic_location -= rtx.coding_exons()
        regulatory_location = intragenic_location

        # Upstream region ...
        if upstream_extension > 0:
            upstream_location = tx.tss_shifted_1bp_upstream().extend_upstream(
                upstream_extension, chromosome2length
            )
            # Changed on 02/09/2011: check if upstream location is not empty ...
            if not upstream_location.isempty():
                for rtx in Transcript.load_by_location(
                    connection,
                    upstream_location.chromosome,
                    upstream_location.span(),
                    only_NMs=is_NM,
                ):
                    upstream_location -= (
                        rtx.coding_exons()
                        if tx.tss_as_bp_location() in rtx.transcript()
                        else rtx.transcript()
                    )
                upstream_location = upstream_location.filter(
                    tx.tss_shifted_1bp_upstream_as_bp_location()
                )
                regulatory_location += upstream_location

        # Downstream region ...
        if downstream_extension > 0:
            downstream_location = tx.tes_shifted_1bp_downstream().extend_downstream(
                downstream_extension, chromosome2length
            )
            # Changed on 02/09/2011: check if downstream location is not empty ...
            if not downstream_location.isempty():
                for rtx in Transcript.load_by_location(
                    connection,
                    downstream_location.chromosome,
                    downstream_location.span(),
                    only_NMs=is_NM,
                ):
                    downstream_location -= (
                        rtx.coding_exons()
                        if tx.tes_as_bp_location() in rtx.transcript()
                        else rtx.transcript()
                    )
                downstream_location = downstream_location.filter(
                    tx.tes_shifted_1bp_downstream_as_bp_location()
                )
                regulatory_location += downstream_location

        gene_name2locations.setdefault(tx.gene_name, []).append(regulatory_location)
    for gene_name in sorted(gene_name2locations.keys()):
        locations = gene_name2locations[gene_name]
        if (
            len({loc.chromosome for loc in locations}) != 1
            or len({loc.on_positive_strand for loc in locations}) != 1
        ):
            print(
                f"Skipped {gene_name:s}: cannot combine regulatory regions.",
                file=sys.stderr,
            )
            continue

        combined_location = None
        for location in locations:
            if combined_location is None:
                combined_location = location
            else:
                combined_location += location

        if combined_location.isempty():
            print(
                f"Skipped {gene_name:s}: no regulatory region remains.",
                file=sys.stderr,
            )
        else:
            for idx, interval in enumerate(combined_location):
                yield (
                    combined_location.chromosome,
                    interval.start,
                    interval.end,
                    f"{gene_name:s}#{idx + 1:d}",
                    "+" if combined_location.on_positive_strand else "-",
                )


def main():
    parser = argparse.ArgumentParser(
        description="""
            Create gene regulatory regions BED file.
        """,
    )

    parser.add_argument(
        "--db",
        dest="database_filename",
        action="store",
        type=str,
        required=True,
        help='Gene SQLite3 database filename, created with "create_genes_database.sh".',
    )
    parser.add_argument(
        "--chrom",
        dest="chromosome_sizes_filename",
        action="store",
        type=str,
        required=True,
        help="Chromosome sizes filename with list of chromosome names in first column (UCSC chrom sizes, FAIDX index or plain list).",
    )
    parser.add_argument(
        "--up",
        dest="upstream_extend",
        action="store",
        type=int,
        required=True,
        help="Extend regulatory region x bp upstream of TSS.",
    )
    parser.add_argument(
        "--down",
        dest="downstream_extend",
        action="store",
        type=int,
        required=True,
        help="Extend regulatory region x bp downstream of TSS.",
    )
    parser.add_argument(
        "--intronic",
        dest="intronic_extend",
        action="store",
        type=int,
        required=True,
        help="Extend regulatory region x bp from introns.",
    )
    parser.add_argument(
        "--mode",
        dest="mode",
        action="store",
        type=str,
        required=True,
        choices={"FullTx", "AllIntrons", "5utrIntron1", "NoTx", "5utr"},
        help="""
            Extend regulatory regions (determined by up, down and intronic arguments) for each gene by including:
              - FullTx: Full transcript (5UTR, all exons and all introns).
              - AllIntrons: All introns (no exons).
              - 5utrIntron1: 5'UTR and 1st intron in CDS
              - NoTx: No transcript.
              - 5utr: 5'UTR.
        """,
    )

    args = parser.parse_args()

    if not os.path.isfile(args.database_filename):
        print(
            """Database file "{args.database_filename}" doesn't exist.""",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.isfile(args.chromosome_sizes_filename):
        print(
            """Chromosome sizes file "{args.chromosome_sizes_filename}" doesn't exist.""",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.mode == "FullTx":
        mode = Modes.FULL_TRANSCRIPT
    elif args.mode == "AllIntrons":
        mode = Modes.ALL_INTRONS
    elif args.mode == "5utrIntron1":
        mode = Modes.UTR5_INTRON1
    elif args.mode == "NoTx":
        mode = Modes.NO_TRANSCRIPT
    elif args.mode == "5utr":
        mode = Modes.UTR5

    chromosomes = set(
        chromosome_iterator(
            filename=args.chromosome_sizes_filename,
            exclude_alt=True,
            exclude_random=True,
        )
    )
    with sqlite3.connect(args.database_filename) as connection:
        chromosome2length = load_chromosome_lengths(connection)
        for columns in regulatory_regions_iterator(
            connection,
            chromosomes,
            chromosome2length,
            args.upstream_extend,
            args.downstream_extend,
            args.intronic_extend,
            mode,
        ):
            print("{0:s}\t{1:d}\t{2:d}\t{3:s}\t0\t{4:s}".format(*columns))


if __name__ == "__main__":
    main()

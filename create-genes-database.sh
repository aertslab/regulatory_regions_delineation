#!/usr/bin/env bash

if [ ${#@} -ne 3 ]; then
	printf 'Wrong number of input arguments.\n\n';
	printf 'Usage:\n    bash create-genes-database.sh <refgene_table_filename> <genome_2bit_filename> <database_name>\n\n';
	printf 'Link to input files:\n';
	printf '    - refgene_table_filename:  ftp://hgdownload.cse.ucsc.edu/goldenPath/<assembly>/database/refGene.txt.gz\n';
	printf '    - genome_2bit_filename:    ftp://hgdownload.cse.ucsc.edu/goldenPath/<assembly>/bigZips/<assembly>.2bit\n\n';
	exit 2;
fi

REFGENE_TABLE_FILENAME="${1}";
GENOME_2BIT_FILENAME="${2}";
DATABASE_FILENAME="${3}";


# Check if all input files exist and have the right extension,
# so there might be a high chance that the correct input files are used.

if [ ! -f "${REFGENE_TABLE_FILENAME}" ] ; then
    printf 'ERROR: RefGene Table file "%s" could not be found.\n' "${REFGENE_TABLE_FILENAME}";
    exit 1;
fi

if [ "${REFGENE_TABLE_FILENAME}" = "${REFGENE_TABLE_FILENAME%.txt.gz}" ] ; then
    printf 'ERROR: RefGene Table file "%s" should end with ".txt.gz".\n' "${REFGENE_TABLE_FILENAME}";
    exit 1;
fi

if [ ! -f "${GENOME_2BIT_FILENAME}" ] ; then
    printf 'ERROR: Genome 2bit file "%s" could not be found.\n' "${GENOME_2BIT_FILENAME}";
    exit 1;
fi

if [ "${GENOME_2BIT_FILENAME}" = "${GENOME_2BIT_FILENAME%.2bit}" ] ; then
    printf 'ERROR: Genome 2bit file "%s" should end with ".2bit".\n' "${GENOME_2BIT_FILENAME}";
    exit 1;
fi

if [ $(type sqlite3 > /dev/null 2>&1; echo $?;) -ne 0 ] ; then
    printf 'ERROR: Add "sqlite3" to ${PATH}.\n';
    exit 1;
fi

if [ $(type twoBitInfo > /dev/null 2>&1; echo $?;) -ne 0 ] ; then
    printf 'ERROR: Add "twoBitInfo" to ${PATH}.\n';
    exit 1;
fi


# Specify all SQL queries.

# Original SQL scheme for refGene table can be found at:
#   ftp://hgdownload.cse.ucsc.edu/goldenPath/<assembly>/database/refGene.sql

CREATE_TABLE_GENES_QUERY='
CREATE TABLE `genes` (
    `geneID` VARCHAR(255) NOT NULL,
    `chromosome` VARCHAR(255) NOT NULL,
    `strand` CHAR(1) NOT NULL,
    `txStart` INTEGER NOT NULL,
    `txEnd` INTEGER NOT NULL,
    `cdsStart` INTEGER NOT NULL,
    `cdsEnd` INTEGER NOT NULL,
    `exonCount` INTEGER NOT NULL,
    `exonStarts` LONGBLOB NOT NULL,
    `exonEnds` LONGBLOB NOT NULL,
    `geneName` VARCHAR(255) NOT NULL
);
'

CREATE_INDICES_FOR_GENES_TABLE_QUERY='
CREATE INDEX `geneID` ON `genes` (`geneID`);
CREATE INDEX `geneName` ON `genes` (`geneName`);
CREATE INDEX `location` ON `genes` (`chromosome`, `strand`, `txStart`, `txEnd`);
';

CREATE_TABLE_CHROMOSOMES_QUERY='
CREATE TABLE `chromosomes` (
    `name` VARCHAR(255) NOT NULL,
    `length` INTEGER NOT NULL
);'

CREATE_INDEX_FOR_CHROMOSOMES_TABLE_QUERY='
CREATE INDEX `name` ON `chromosomes` (`name`);
'


# Remove old database file, if it exists.
rm -f "${DATABASE_FILENAME}";


# Create "genes" and "chromsomes" tables.
sqlite3 -line \
    "${DATABASE_FILENAME}" \
    "${CREATE_TABLE_GENES_QUERY}${CREATE_TABLE_CHROMOSOMES_QUERY}"


# Import refGene data in "genes" table:
#   - Only keep the following columns:
#       name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, name2
zcat "${REFGENE_TABLE_FILENAME}" \
    | awk -F '\t' -v 'OFS=\t' '{ print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $13; }' \
    | sqlite3 \
        -separator $'\t' \
        -line \
        "${DATABASE_FILENAME}" \
        '.import /dev/stdin genes';


# Import chromosome name and chromosome size data in "chromosomes" tables.
twoBitInfo "${GENOME_2BIT_FILENAME}" stdout \
    | sqlite3 \
        -separator $'\t' \
        -line \
        "${DATABASE_FILENAME}" \
        '.import /dev/stdin chromosomes';


# Create indices for "genes" and "chromosomes" tables.
sqlite3 \
    -line \
    "${DATABASE_FILENAME}" \
    "${CREATE_INDICES_FOR_GENES_TABLE_QUERY}${CREATE_INDEX_FOR_CHROMOSOMES_TABLE_QUERY}";

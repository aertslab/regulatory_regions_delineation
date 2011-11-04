#!/usr/bin/env bash

if [ $# -ne 3 ]; then
	echo "Wrong number of input arguments."
	echo
	echo "Usage: bash create-genes-database.sh <table_filename> <genome_2bit_filename> <database_name>"
	exit 2
fi
TABLE_FILENAME=$1
GENOME_FILENAME=$2
DATABASE_FILENAME=$3

rm -f ${DATABASE_FILENAME}
# Only keep columns: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, name2
cat ${TABLE_FILENAME} | sed '1d' | cut -f2,3,4,5,6,7,8,9,10,11,13 | tr '\t' '|' > sqlite3-genes-import.psv
# Find length of chromosomes ...
twoBitInfo ${GENOME_FILENAME} stdout | tr '\t' '|' > sqlite3-chromosomes-import.psv
sqlite3 -line ${DATABASE_FILENAME} "create table genes ( ID varchar(255), \
chromosome varchar(255), strand char(1), txStart integer, txEnd integer, \
cdsStart integer, cdsEnd integer, exonCount integer, exonStarts blob, exonEnds blob, altName varchar(255));"
sqlite3 -separator '|' -line ${DATABASE_FILENAME} ".import sqlite3-genes-import.psv genes"
sqlite3 -line ${DATABASE_FILENAME} "create index ID on genes ( ID );"
sqlite3 -line ${DATABASE_FILENAME} "create index altName on genes ( altName );"
sqlite3 -line ${DATABASE_FILENAME} "create index location on genes ( chromosome, strand, txStart, txEnd );"
sqlite3 -line ${DATABASE_FILENAME} "create table chromosomes ( name varchar(255), length integer);"
sqlite3 -separator '|' -line ${DATABASE_FILENAME} ".import sqlite3-chromosomes-import.psv chromosomes"
sqlite3 -line ${DATABASE_FILENAME} "create index name on chromosomes ( name );"
rm -f sqlite3-genes-import.psv sqlite3-chromosomes-import.psv
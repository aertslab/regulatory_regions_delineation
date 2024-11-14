#!/usr/bin/env bash


create_genes_database () {
    if [ ${#@} -ne 3 ]; then
        printf 'Wrong number of input arguments.\n\n';
        printf 'Usage:\n    bash create_genes_database.sh <refgene_table_filename|gene_pred_filename> <genome_2bit_filename> <gene_sqlite3_database_filename>\n\n';
        printf 'Arguments:\n';
        printf '  - refgene_table_filename:\n';
        printf '      Download premade refGene files from:\n';
        printf '        ftp://hgdownload.cse.ucsc.edu/goldenPath/<assembly>/database/refGene.txt.gz\n';
        printf '  - gene_pred_filename:\n';
        printf '      Create Gene Predictions (Extended) file from GTF/GFF with:\n';
        printf '        - gtfToGenePred -genePredExt genes.gtf /dev/stdout | gzip -c > genes.txt.gz\n';
        printf '          (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred)\n';
        printf '        - gtf2tbl.py genes.gtf | gzip -c > genes.txt.gz\n';
        printf '        - gff2tbl.py genes.gff | gzip -c > genes.txt.gz\n';
        printf '  - genome_2bit_filename:\n';
        printf '      Download genome 2bit files from:\n';
        printf '        ftp://hgdownload.cse.ucsc.edu/goldenPath/<assembly>/bigZips/<assembly>.2bit\n';
        printf '      Create genome 2bit files from FASTA files with:\n';
        printf '        faToTwoBit genome.fa genome.2bit \n';
        printf '        (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit)\n';
        printf '  - gene_sqlite3_database_filename:\n';
        printf '      Output gene SQLite3 database filename.\n\n';
        return 2;
    fi

    local gene_pred_filename="${1}";
    local genome_2bit_filename="${2}";
    local gene_sqlite3_database_filename="${3}";


    # Check if all input files exist and have the right extension,
    # so there might be a high chance that the correct input files are used.

    if [ ! -f "${gene_pred_filename}" ] ; then
        printf 'ERROR: GenePred file or RefGene Table file "%s" could not be found.\n' "${gene_pred_filename}";
        return 1;
    fi

    if [ "${gene_pred_filename}" = "${gene_pred_filename%.txt.gz}" ] ; then
        printf 'ERROR: GenePred file or RefGene Table file "%s" should end with ".txt.gz".\n' "${gene_pred_filename}";
        return 1;
    fi

    if [ ! -f "${genome_2bit_filename}" ] ; then
        printf 'ERROR: Genome 2bit file "%s" could not be found.\n' "${genome_2bit_filename}";
        return 1;
    fi

    if [ "${genome_2bit_filename}" = "${genome_2bit_filename%.2bit}" ] ; then
        printf 'ERROR: Genome 2bit file "%s" should end with ".2bit".\n' "${genome_2bit_filename}";
        return 1;
    fi

    if [ $(type sqlite3 > /dev/null 2>&1; echo $?;) -ne 0 ] ; then
        printf 'ERROR: Add "sqlite3" to ${PATH}.\n';
        return 1;
    fi

    if [ $(type twoBitInfo > /dev/null 2>&1; echo $?;) -ne 0 ] ; then
        printf 'ERROR: Add "twoBitInfo" (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo) to ${PATH}.\n';
        return 1;
    fi


    # Specify all SQL queries.

    # Original SQL scheme for refGene table can be found at:
    #   ftp://hgdownload.cse.ucsc.edu/goldenPath/<assembly>/database/refGene.sql

    local CREATE_TABLE_GENES_QUERY='
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

    local CREATE_INDICES_FOR_GENES_TABLE_QUERY='
    CREATE INDEX `geneID` ON `genes` (`geneID`);
    CREATE INDEX `geneName` ON `genes` (`geneName`);
    CREATE INDEX `location` ON `genes` (`chromosome`, `strand`, `txStart`, `txEnd`);
    ';

    local CREATE_TABLE_CHROMOSOMES_QUERY='
    CREATE TABLE `chromosomes` (
        `name` VARCHAR(255) NOT NULL,
        `length` INTEGER NOT NULL
    );'

    local CREATE_INDEX_FOR_CHROMOSOMES_TABLE_QUERY='
    CREATE INDEX `name` ON `chromosomes` (`name`);
    '


    # Remove old database file, if it exists.
    rm -f "${gene_sqlite3_database_filename}";


    # Create "genes" and "chromsomes" tables.
    sqlite3 -line \
        "${gene_sqlite3_database_filename}" \
        "${CREATE_TABLE_GENES_QUERY}${CREATE_TABLE_CHROMOSOMES_QUERY}"


    # Import refGene data in "genes" table:
    #   - Only keep the following columns:
    #       name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, name2
    zcat "${gene_pred_filename}" \
        | awk \
            -F '\t' \
            -v 'OFS=\t' \
            '
            {
                if ($1 ~ /^#/) {
                    # Skip header line.
                    next;
                } else if (NF == 15) {
                    # GenePred file without "bin" column.
                    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $12;
                } else if (NF == 16) {
                    # GenePred file with "bin" column.
                    print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $13;
                }
            }' \
        | sqlite3 \
            -separator $'\t' \
            -line \
            "${gene_sqlite3_database_filename}" \
            '.import /dev/stdin genes';


    # Import chromosome name and chromosome size data in "chromosomes" tables.
    twoBitInfo "${genome_2bit_filename}" stdout \
        | sqlite3 \
            -separator $'\t' \
            -line \
            "${gene_sqlite3_database_filename}" \
            '.import /dev/stdin chromosomes';


    # Create indices for "genes" and "chromosomes" tables.
    sqlite3 \
        -line \
        "${gene_sqlite3_database_filename}" \
        "${CREATE_INDICES_FOR_GENES_TABLE_QUERY}${CREATE_INDEX_FOR_CHROMOSOMES_TABLE_QUERY}";
}



create_genes_database "${@}";

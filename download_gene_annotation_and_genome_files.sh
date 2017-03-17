#!/usr/bin/env bash

ICISTARGET_DATA_DIR='/media/data/lcb/icistarget/data';

# Default UCSC download server.
#UCSC_DOWNLOAD_MIRROR='hgdownload.cse.ucsc.edu';
# San Diego UCSC download server.
UCSC_DOWNLOAD_MIRROR='hgdownload-sd.cse.ucsc.edu';


download_refseq_gene_annotation_file () {
    local refseq_version="${1}";
    local species="${2}";
    local assembly="${3}";
    local refseq_basename="${4}";
    refseq_basename="${refseq_basename:=refGene}";

    printf "\n\nDownloading refseq gene annotation to:\n  %s\nfor %s ...\n\n" \
        "${ICISTARGET_DATA_DIR}/annotations/${species}/${assembly}/" \
        "${species} ${assembly} ${refseq_basename}";

    mkdir -p "${ICISTARGET_DATA_DIR}/annotations/${species}/${assembly}/";

    # RefSeq schema.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/database/${refseq_basename}.sql" \
        "${ICISTARGET_DATA_DIR}/annotations/${species}/${assembly}/${refseq_basename}.${assembly}__refseq-r${refseq_version}.sql";

    # RefSeq annotation.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/database/${refseq_basename}.txt.gz" \
        "${ICISTARGET_DATA_DIR}/annotations/${species}/${assembly}/${refseq_basename}.${assembly}__refseq-r${refseq_version}.txt.gz";
}


download_2bit_file () {
    local assembly="${1}";

    printf "\n\nDownloading 2bit file to:\n  %s\nfor %s ...\n\n" \
        "${ICISTARGET_DATA_DIR}/genomes/2bit/" \
        "${assembly}";

    mkdir -p "${ICISTARGET_DATA_DIR}/genomes/2bit/";

    # Download 2bit file for specified assembly.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/bigZips/${assembly}.2bit" \
        "${ICISTARGET_DATA_DIR}/genomes/2bit/";

    # Download chrom sizes file for specified assembly.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/bigZips/${assembly}.chrom.sizes" \
        "${ICISTARGET_DATA_DIR}/genomes/chrom_sizes/";
}


download_chain_file () {
    local from_assembly="${1}";
    local to_assembly="${2}";

    printf "\n\nDownloading chain file to:\n  %s\nfor %s to %s ...\n\n" \
        "${ICISTARGET_DATA_DIR}/genomes/liftOver/" \
        "${from_assembly}" \
        "${to_assembly}";

    mkdir -p "${ICISTARGET_DATA_DIR}/genomes/liftOver/";

    # Download chain file for specified from_assembly and to_assembly.
    # (Convert first character of to_assembly variable to uppercase with ${to_assembly^}.)
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${from_assembly}/liftOver/${from_assembly}To${to_assembly^}.over.chain.gz" \
        "${ICISTARGET_DATA_DIR}/genomes/liftOver/";
}


download_2bit_files_and_chain_files_for_liftover_genomes_filename () {
    # Liftover genomes filename:
    #   - column 1: to_assembly.
    #   - column 2: full species name of to_assembly.
    #   - column 3: from_assembly (should be the same for the whole file).
    local liftover_genomes_filename="${1}";

    local assembly='';
    local to_assembly='';

    # Get "from assembly" from the third column.
    local from_assembly=$(cut -f3 "${liftover_genomes_filename}" | sort -u);

    # Get "to assemblies" from the first column.
    local to_assemblies=$(cut -f1 "${liftover_genomes_filename}");

    # Download all needed 2bit files.
    for assembly in ${from_assembly} ${to_assemblies}; do
        download_2bit_file "${assembly}";
    done

    # Download all needed chain files.
    for to_assembly in ${to_assemblies}; do
        download_chain_file "${from_assembly}" "${to_assembly}";
    done
}


main () {
    # Download refSeq version 80 gene annotation files for the following species.
    download_refseq_gene_annotation_file 80 'danio_rerio' 'danRer10'
    download_refseq_gene_annotation_file 80 'mus_musculus' 'mm10'
    download_refseq_gene_annotation_file 80 'homo_sapiens' 'hg38' 'ncbiRefSeqCurated'

    # Download all needed 2bit and chain files.
    download_2bit_files_and_chain_files_for_liftover_genomes_filename "configurations/danRer10.liftover_genomes.lst";
    download_2bit_files_and_chain_files_for_liftover_genomes_filename "configurations/hg38.liftover_genomes.lst";
    download_2bit_files_and_chain_files_for_liftover_genomes_filename "configurations/mm10.liftover_genomes.lst";
}

#main;

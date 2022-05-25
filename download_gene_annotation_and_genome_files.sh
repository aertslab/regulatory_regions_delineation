#!/usr/bin/env bash

DATA_DIR='./data';

# Default UCSC download server.
UCSC_DOWNLOAD_MIRROR='hgdownload.cse.ucsc.edu';


get_download_command () {
    if ( ! [ -z "${download_command}" ] && type "${download_command}" > /dev/null 2>&1 ) ; then
        true;
    elif type curl > /dev/null 2>&1 ; then
        download_command='curl';
    elif type wget > /dev/null 2>&1 ; then
        download_command='wget';
    else
        printf 'Error: "cURL" and "wget" are not found in the path.\n\n' > 2;
        printf '       Install one of them or specify the full path to the binary:\n\n' > 2;
        printf '           download_command=/my/custom/path/to/curl\n\n' > 2;
        return 1;
    fi

    printf "Using \"%s\" to download files: download_command='%s'\n" \
        "${download_command}" \
        "${download_command}" > 2;
}


get_download () {
    local download_url="${1}";
    local download_filename="${2}";
    local hide_download_progress="${3}";

    if [ ${#@} -lt 2 ] ; then
        printf 'Usage: get_download download_url download_filename hide_download_progress=[01]\n\n';
        printf 'Examples:\n';
        printf '  Print content of url to standard output and print download progress:\n';
        printf '    download_to_stdout http://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/ stdout\n\n';
        printf '  Write content of url to a filename and print download progress:\n';
        printf '    download_to_stdout http://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/ release-notes-index.txt\n\n';
        printf '  Print content of url to stdout and do not print download progress:\n';
        printf '    download_to_stdout http://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/ stdout 1\n\n';
        return 1;
    fi

    if [ "${download_filename}" = "stdout" ] ; then
        download_filename='-';
    fi

    if [ -z "${download_command}" ] ; then
        # Check if we find cURL of wget (else exit).
        get_download_command || return 1;
    fi

    case "${download_command}" in
        *curl)
            # Download with cURL.
            if [ "${hide_download_progress}" = 1 ] ; then
                "${download_command}" \
                    -R \
                    -L \
                    -o "${download_filename}" \
                    "${download_url}" 2> /dev/null;
            else
                "${download_command}" \
                    -R \
                    -L \
                    -o "${download_filename}" \
                    "${download_url}";
            fi
            ;;
        *wget)
            # Download with wget.
            if [ "${hide_download_progress}" = 1 ] ; then
                "${download_command}" \
                    -O "${download_filename}" \
                    "${download_url}" 2> /dev/null;
            else
                "${download_command}" \
                    -O "${download_filename}" \
                    "${download_url}";
            fi
            ;;
        *)
            printf 'Error: Unknown download command "%s". Expected "curl" or "wget".\n' \
                "${download_command}" > 2;
            return 1;
            ;;
    esac;
}


get_last_refseq_version () {
    # Get latest released version of refSeq by looking at "RefSeq-release-*.txt"
    # in http://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/.
    # In general UCSC keeps gene annotation up to date with refSeq releases.
    refseq_version=$(
        get_download http://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/ stdout 1 \
            | sed -n 's/.\+RefSeq-release\([0-9]\+\).\+/\1/p'
    );

    echo "${refseq_version}";
}


download_refseq_gene_annotation_file () {
    if [ ${#@} -ne 3 ] ; then
        printf 'Usage: download_refseq_gene_annotation_file refseq_version assembly <refGene|ncbiRefSeq>\n\n';
        printf 'Parameters:\n';
        printf '  - refseq_version:\n';
        printf '      Specify refSeq version used on UCSC: e.g.: 90\n';
        printf '      If refSeq version is set to 0, the last refSeq version will be used.\n';
        printf '      UCSC normally updates gene annotation to the last version of refSeq,\n';
        printf '      so in general you should just set this parameter to 0.\n\n';
        printf '  - assembly:\n';
        printf '      Specify assembly version: e.g. hg38\n';
        printf '      Check the following links for possible assembly version names:\n';
        printf '        - http://hgdownload.cse.ucsc.edu/downloads.html\n';
        printf '        - http://hgdownload.cse.ucsc.edu/goldenPath/\n\n';
        printf '  - refGene|ncbiRefSeq:\n';
        printf '      Download refGene or ncbiRefSeq gene annotation.\n\n';
        printf 'Examples:\n';
        printf '  Download last ncbiRefSeqCurated gene annotation for Homo sapiens (hg38) from UCSC:\n\n';
        printf '  download_refseq_gene_annotation_file 0 hg38 ncbiRefSeqCurated\n\n';
        return 1;
    fi

    local -i refseq_version="${1}";
    local assembly="${2}";
    local refseq_basename="${3}";
    refseq_basename="${refseq_basename:=refGene}";

    if [ "${refseq_version}" -eq 0 ] ; then
        printf 'Getting last refSeq version ...\n';
        refseq_version=$(get_last_refseq_version);
        printf 'Setting refSeq version: %d\n' "${refseq_version}";
    fi

    printf "\n\nDownloading refseq gene annotation to:\n  %s\nfor %s ...\n\n" \
        "${DATA_DIR}/annotations/${assembly}/" "${assembly}";

    mkdir -p "${DATA_DIR}/annotations/${assembly}/";

    # RefSeq schema.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/database/${refseq_basename}.sql" \
        "${DATA_DIR}/annotations/${assembly}/${refseq_basename}.${assembly}__refseq-r${refseq_version}.sql";

    # RefSeq annotation.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/database/${refseq_basename}.txt.gz" \
        "${DATA_DIR}/annotations/${assembly}/${refseq_basename}.${assembly}__refseq-r${refseq_version}.txt.gz";
}


download_2bit_file () {
    if [ ${#@} -ne 1 ] ; then
        printf 'Usage: download_2bit_file assembly\n';
        return 1;
    fi

    local assembly="${1}";

    printf "\n\nDownloading 2bit file to:\n  %s\nfor %s ...\n\n" \
        "${DATA_DIR}/genomes/2bit/" \
        "${assembly}";

    mkdir -p "${DATA_DIR}/genomes/2bit/";

    # Download 2bit file for specified assembly.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/bigZips/${assembly}.2bit" \
        "${DATA_DIR}/genomes/2bit/";

    # Download chrom sizes file for specified assembly.
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${assembly}/bigZips/${assembly}.chrom.sizes" \
        "${DATA_DIR}/genomes/chrom_sizes/";
}


download_chain_file () {
    if [ ${#@} -ne 2 ] ; then
        printf 'Usage: download_chain_file from_assembly to_assembly\n';
        return 1;
    fi

    local from_assembly="${1}";
    local to_assembly="${2}";

    printf "\n\nDownloading chain file to:\n  %s\nfor %s to %s ...\n\n" \
        "${DATA_DIR}/genomes/liftOver/" \
        "${from_assembly}" \
        "${to_assembly}";

    mkdir -p "${DATA_DIR}/genomes/liftOver/";

    # Download chain file for specified from_assembly and to_assembly.
    # (Convert first character of to_assembly variable to uppercase with ${to_assembly^}.)
    rsync -avzP \
        "rsync://${UCSC_DOWNLOAD_MIRROR}/goldenPath/${from_assembly}/liftOver/${from_assembly}To${to_assembly^}.over.chain.gz" \
        "${DATA_DIR}/genomes/liftOver/";
}


download_2bit_files_and_chain_files_for_liftover_genomes_filename () {
    if [ ${#@} -ne 1 ] ; then
        printf 'Usage: download_2bit_files_and_chain_files_for_liftover_genomes_filename liftover_genomes_filename\n\n';
        printf 'Liftover genomes filename: TAB-separated file with 3 columns\n';
        printf '  - column 1: to_assembly: assembly name to which you want to convert.\n';
        printf '  - column 2: full species name of to_assembly.\n';
        printf '  - column 3: from_assembly (should be the same for the whole file):\n';
        printf '              assembly name from which you want to convert.\n\n';
        return 1;
    fi

    # Liftover genomes filename:
    #   - column 1: to_assembly.
    #   - column 2: full species name of to_assembly.
    #   - column 3: from_assembly (should be the same for the whole file).
    local liftover_genomes_filename="${1}";

    if [ ! -f "${liftover_genomes_filename}" ] ; then
        printf 'Error: LiftOver genomes filename "%s" could not be found.\n' \
            "${liftover_genomes_filename}";
        return 1;
    fi

    local assembly='';
    local to_assembly='';

    # Get "from assembly" from the third column.
    local from_assembly=$(cut -f3 "${liftover_genomes_filename}" | sort -u);

    # Get "to assemblies" from the first column.
    local to_assemblies=$(cut -f1 "${liftover_genomes_filename}");

    # Download all needed 2bit files.
    for assembly in ${from_assembly} ${to_assemblies} ; do
        download_2bit_file "${assembly}";
    done

    # Download all needed chain files.
    for to_assembly in ${to_assemblies} ; do
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

show_help () {
    printf 'Import functions in this script:\n\n';
    printf '  source %s\n\n' "${BASH_SOURCE[0]}";

    printf 'Specify DATA_DIR to which files will be downloaded:\n\n';
    printf '  DATA_DIR=%s\n\n' "${DATA_DIR}";

    printf 'Functions:\n';
    printf '  - show_help: Show this help text.\n';
    printf '  - get_download_command: Find cURL or wget.\n';
    printf '  - get_download: Download a file from an URL.\n';
    printf '  - get_last_refseq_version: Get last refSeq version.\n';
    printf '  - download_refseq_gene_annotation_file: Download refSeq gene annotation file.\n';
    printf '  - download_2bit_files_and_chain_files_for_liftover_genomes_filename: Download 2bit and chain files.\n\n';
}


show_help;

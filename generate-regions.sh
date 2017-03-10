#!/usr/bin/env bash

# Dependencies:
#	create-genes-database.sh, generate-liftover-fasta.sh
#   create-regulatory-regions-bed.py, interval.py, piecewiselocation.py, transcript.py
#	twoBitInfo, twoBitToFa, sqlite3


# Parameters ...
if [ $# -ne 2 ]; then
	printf 'Usage: ./generate-regions.sh ini_file output_dir\n\nERROR: Wrong number of arguments.\n\n' > /dev/stderr;
	exit 2;
fi

ini_filename="${1}";
output_dir="${2}";

if [ ! -f "${ini_filename}" ] ; then
	printf 'ERROR: Ini file does not exist.\n';
	exit 2
fi

if [ ! -d "${output_dir}" ] ; then
	mkdir "${output_dir}";
fi


get_config_parameter () {
    local ini_filename="${1}";
    local config_parameter="${2}";

    awk -F '=' -v config_parameter="${config_parameter}" \
        '{ if ($1 == config_parameter) { print $2; } }' \
        "${ini_filename}";
}


file_does_not_exists () {
    local filename="${1}";

    if [ ! -f "${filename}" ] ; then
        printf '\nERROR: File "%s" does not exist.\n\n' "${filename}" > /dev/stderr;
        exit 1;
    fi
}


# Parse config file.

# Get i-cisTarget data directory which will be prepended to various file paths.
icistarget_data_dir=$(get_config_parameter "${ini_filename}" 'icistarget_data_dir');

ucsc_table="${icistarget_data_dir}"/$(get_config_parameter "${ini_filename}" 'ucsc_table');
id_description_table="${icistarget_data_dir}"/$(get_config_parameter "${ini_filename}" 'id_description_table');
base_genome_2bit_file="${icistarget_data_dir}"/$(get_config_parameter "${ini_filename}" 'base_genome_2bit_file');
genomes_liftover_table="${icistarget_data_dir}"/$(get_config_parameter "${ini_filename}" 'genomes_liftover_table');

chromosomes=$(get_config_parameter "${ini_filename}" 'chromosomes');
base_genome_id=$(get_config_parameter "${ini_filename}" 'base_genome_id');
delineation=$(get_config_parameter "${ini_filename}" 'delineation');

# Get upstream, intronic and downstream extension from config file or initialize to 0 if not found.
upstream_extension_in_bp=$(get_config_parameter "${ini_filename}" 'upstream_extension_in_bp');
upstream_extension_in_bp="${upstream_extension_in_bp:=0}";
intronic_extension_in_bp=$(get_config_parameter "${ini_filename}" 'downstream_extension_in_bp');
intronic_extension_in_bp="${intronic_extension_in_bp:=0}";
downstream_extension_in_bp=$(get_config_parameter "${ini_filename}" 'downstream_extension_in_bp');
downstream_extension_in_bp="${downstream_extension_in_bp:=0}";

# Get number of species which will be used for doing liftOver.
nbr_of_species=$(
    awk -F '\t' -v base_genome_id="${base_genome_id}" \
        '{ if ($3 == base_genome_id) { nbr_of_species += 1 } } END { print nbr_of_species }' \
        "${genomes_liftover_table}"
);


echo "Configuration settings";
echo "----------------------";
echo;
echo "UCSC gene table              = ${ucsc_table}";
echo "Gene ID to description table = ${id_description_table}";
echo "Base genome 2bit file        = ${base_genome_2bit_file}";
echo "Base genome ID               = ${base_genome_id}";
echo "Liftover genomes table       = ${genomes_liftover_table}";
echo "Chromosomes                  = ${chromosomes}";
echo "Upstream extend              = ${upstream_extension_in_bp}bp";
echo "Delineation                  = ${delineation}";
echo "Intronic extend              = ${intronic_extension_in_bp}bp";
echo "Downstream extend            = ${downstream_extension_in_bp}bp";
echo "Included (${nbr_of_species}) species:"
awk -F '\t' -v base_genome_id="${base_genome_id}" \
    '{ if ($3 == base_genome_id) { print "\t" $1 "\t" $2; } }' \
    "${genomes_liftover_table}";
echo;


# Check if some of the necessary files exist.
file_does_not_exists "${ucsc_table}";
file_does_not_exists "${id_description_table}";
file_does_not_exists "${base_genome_2bit_file}";
file_does_not_exists "${genomes_liftover_table}";


create_region_description() {
	local upstream_extension_in_bp="${1}";
	local delineation="${2}";
	local downstream_extension_in_bp="${3}";
	local intronic_extension_in_bp="${4}";

    case "${delineation}" in
        'FullTx')
            suffix='-full-transcript';;
        'AllIntrons')
            suffix='-introns';;
        'NoTx')
            suffix='';;
        '5utr')
            suffix='-5utr';;
        '5utr-intron1')
            suffix='-5utr-intron1';;
        *)
            printf '\nERROR: Invalid delineation "%s".\n\n' "${delineation}" > /dev/stderr;
	        exit 2;;
    esac

	if [ ${intronic_extension_in_bp} -gt 0 ] ; then
	    intronic="-tss-downstream${intronic_extension_in_bp}";
	else
	    intronic='';
	fi

	if [ ${downstream_extension_in_bp} -gt 0 ] ; then
		downstream="-downstream${downstream_extension_in_bp}";
	else
	    downstream="";
	fi

	if [ ${upstream_extension_in_bp} -gt 0 ] ; then
		upstream="-upstream${upstream_extension_in_bp}";
	else
	    upstream="";
	fi

	if [ ${downstream_extension_in_bp} -gt 0 -o ${upstream_extension_in_bp} -gt 0 ]; then
		echo "-limited${upstream}${downstream}${intronic}${suffix}";
	else
	    echo "${intronic}${suffix}";
	fi
}


create_filename() {
	local base_genome_id="${1}";
	local output_dir="${2}";
	local upstream_extension_in_bp="${3}";
	local delineation="${4}";
	local downstream_extension_in_bp="${5}";
	local extension="${6}";

	region_description=$(
	    create_region_description \
	        "${upstream_extension_in_bp}" \
	        "${delineation}" \
	        "${downstream_extension_in_bp}" \
	        "${intronic_extension_in_bp}"
	);

	echo "${output_dir}/${base_genome_id}${region_description}.${extension}";
}


# Create SQLite3 gene database ...
printf '\nCreate SQLite3 gene database ...\n\n';
db_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "sqlite3.db")
bash create-genes-database.sh ${ucsc_table} ${base_genome_2bit_file} ${db_filename}

# Create BED file ...
printf '\nCreate BED file ...\n\n';
bed_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "bed")
echo "${chromosomes}" | tr ';' '\n' > $output_dir/chromosomes.tmp
python create-regulatory-regions-bed.py ${db_filename} $output_dir/chromosomes.tmp ${upstream_extension_in_bp} ${downstream_extension_in_bp} ${intronic_extension_in_bp} ${delineation} > ${bed_filename}
rm -f $output_dir/chromosomes.tmp

# Create distribution of length of search space ...
figure_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "png")
/.generate-histogram.sh ${bed_filename}

# Create gene ID file and gene description table ...
printf '\nCreate gene ID file and gene description table ...\n\n';
id_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "gene-ids")
cat ${bed_filename} | cut -f4 | cut -d'#' -f1 | sort -u > ${id_filename}
description_table_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "gene-descriptions")
cat ${id_description_table} | sort -k1 | join -a 1 -1 1 -2 1 ${id_filename} - | sed -e 's/ /|/' | tr '|' '\t' > ${description_table_filename}
regions_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "regions")
cat ${bed_filename} | cut -f4 | sort -u > ${regions_filename}
region_table_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "gene-regions")
python generate-gene-region-table.py ${bed_filename} > ${region_table_filename}

# Analysis of lost gene IDs ...
printf '\nAnalysis of lost gene IDs ...\n\n';
lostids_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "lost-gene-ids")
cat ${ucsc_table} | sed '1d' | cut -f13 | sort -u > $output_dir/geneids-all.tmp
grep -vxF -f ${id_filename} $output_dir/geneids-all.tmp > ${lostids_filename}
rm -f $output_dir/geneids-all.tmp
nr_of_lost_genes=`cat ${lostids_filename} | wc -l`
echo "Lost ${nr_of_lost_genes}# genes ..."

# Create i-cisTarget ini file ...
printf '\nCreate i-cisTarget ini file ...\n\n';
full_output_dir=`cd ${output_dir}; pwd`
full_regions_filename=`python -c "import os.path; print os.path.realpath('${regions_filename}');"`
full_region_table_filename=`python -c "import os.path; print os.path.realpath('${region_table_filename}');"`
cat cistargetx-template.ini \
		| sed -e "s!HOME!${HOME}!; s!DATADIR!${DATADIR}!; s!OUTPUTDIR!${full_output_dir}!; s!REGIONSFILE!${full_regions_filename}!; s!LUTFILE!${full_region_table_filename}!" \
		> ${output_dir}/cistargetx-template.ini

# Create FASTA file ...
printf '\nCreate FASTA file ...\n\n';
fasta_filename=$(create_filename ${base_genome_id} ${output_dir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "fasta")
twoBitToFa -bed=${bed_filename} ${base_genome_2bit_file} ${fasta_filename}
#Usage of -bed instead of -seqList is easier and safe ...
#awk '{printf "%s:%d-%d\n", $1, $2, $3}' ${bed_filename} | twoBitToFa ${base_genome_2bit_file} ${fasta_filename} -seqList=stdin

# Liftover procedure ...
printf '\nLiftover procedure ...\n\n';
${DATADIR}/genomes/generate-liftover-fasta.sh ${bed_filename} \
		${base_genome_id} ${output_dir} $(create_region_description ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} ${intronic_extension_in_bp}) ${genomes_liftover_table}
		
# Copy ini file to output folder ...
cp ${ini_filename} ${output_dir}

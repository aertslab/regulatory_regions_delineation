#!/usr/bin/env bash

# Dependencies:
#	create-genes-database.sh, generate-liftover-fasta.sh
#   create-regulatory-regions-bed-new.py, interval.py, piecewiselocation.py, transcript.py
#	twoBitInfo, twoBitToFa, sqlite3

# Parameters ...
if [ $# -ne 2 ]; then
	echo "Usage: ./generate-regions.sh ini_file output_dir" > /dev/stderr
	echo "" > /dev/stderr
	echo "Wrong number of arguments." > /dev/stderr
	exit 2
fi
inifile=$1
outputdir=$2

if [ ! -e $inifile ]; then
	echo "Ini file does not exist."
	exit 2
fi

if [ ! -d $outputdir ]; then
	mkdir $outputdir
fi

# Loading ini file ...
ucsc_table=`cat $inifile | grep ucsc_table | cut -d= -f2`
chromosomes=`cat $inifile | grep chromosomes | cut -d= -f2`
id_description_table=`cat $inifile | grep id_description_table | cut -d= -f2`
base_genome_2bit_file=`cat $inifile | grep base_genome_2bit_file | cut -d= -f2`
base_genome_id=`cat $inifile | grep base_genome_id | cut -d= -f2`
genomes_liftover_table=`cat $inifile | grep genomes_liftover_table | cut -d= -f2`
upstream_extension_in_bp=`cat $inifile | grep upstream_extension_in_bp | cut -d= -f2`
downstream_extension_in_bp=`cat $inifile | grep downstream_extension_in_bp | cut -d= -f2`
delineation=`cat $inifile | grep delineation | cut -d= -f2`
nr_of_species=`cat ${genomes_liftover_table} | awk -F '\t' -v r=${base_genome_id} '$3==r {print $1 "\t" $2}' | wc -l`
if [ "$(cat $inifile | grep intronic_extension | wc -l)" -lt "1" ]; then intronic_extension_in_bp=0
else intronic_extension_in_bp=$(cat $inifile | grep intronic_extension | cut -d= -f2); fi

echo "UCSC gene table              = ${ucsc_table}"
echo "Gene ID to description table = ${id_description_table}"
echo "Base genome 2bit file        = ${base_genome_2bit_file}"
echo "Base genome ID               = ${base_genome_id}"
echo "Liftover genomes table       = ${genomes_liftover_table}"
echo "Chromosomes                  = ${chromosomes}"
echo "Upstream extend              = ${upstream_extension_in_bp}bp"
echo "Delineation                  = ${delineation}"
echo "Intronic extend              = ${intronic_extension_in_bp}bp"
echo "Downstream extend            = ${downstream_extension_in_bp}bp"
echo "Included species (${nr_of_species}#): "
cat ${genomes_liftover_table} | awk -F '\t' -v r=${base_genome_id} '$3==r {print $1 "\t" $2}'

create_region_description() {
	local upstream_extension_in_bp=$1
	local delineation=$2
	local downstream_extension_in_bp=$3
	local intronic_extension_in_bp=$4
	
	if [ "${delineation}" == "FullTx" ]; then
		suffix='-full-transcript'
	elif [ "${delineation}" == "AllIntrons" ]; then
		suffix='-introns'
	elif [ "${delineation}" == "NoTx" ]; then
		suffix=''
	elif [ "${delineation}" == "5utr" ]; then
		suffix='-5utr'
	else
		suffix='-5utr-intron1'
	fi

	if [ ${intronic_extension_in_bp} -gt 0 ]; then
	   intronic="-tss-downstream${intronic_extension_in_bp}"
	else intronic=''; fi
	if [ ${downstream_extension_in_bp} -gt 0 ]; then
		downstream="-downstream${downstream_extension_in_bp}"
	else downstream="";	fi
	if [ ${upstream_extension_in_bp} -gt 0 ]; then
		upstream="-upstream${upstream_extension_in_bp}"
	else upstream=""; fi

	if [ ${downstream_extension_in_bp} -gt 0 -o ${upstream_extension_in_bp} -gt 0 ]; then
		echo "-limited${upstream}${downstream}${intronic}${suffix}"
	else echo "${intronic}${suffix}"; fi
}

create_filename() {
	local base_genome_id=$1
	local outputdir=$2
	local upstream_extension_in_bp=$3
	local delineation=$4
	local downstream_extension_in_bp=$5
	local extension=$6
	region_description=$(create_region_description ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} ${intronic_extension_in_bp})
	echo "${outputdir}/${base_genome_id}${region_description}.${extension}"
}

# Create SQLite3 database ...
echo "Create SQLite3 database ..."
db_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "sqlite3.db")
bash create-genes-database.sh ${ucsc_table} ${base_genome_2bit_file} ${db_filename}

# Create BED file ...
echo "Create BED file ..."
bed_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "bed")
echo "${chromosomes}" | tr ';' '\n' > $outputdir/chromosomes.tmp
python create-regulatory-regions-bed-new.py ${db_filename} $outputdir/chromosomes.tmp ${upstream_extension_in_bp} ${downstream_extension_in_bp} ${intronic_extension_in_bp} ${delineation} > ${bed_filename}
rm -f $outputdir/chromosomes.tmp

# Create distribution of length of search space ...
figure_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "png")
/.generate-histogram.sh ${bed_filename}

# Create gene ID file and gene description table ...
echo "Create gene ID file and gene description table ..."
id_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "gene-ids")
cat ${bed_filename} | cut -f4 | cut -d'#' -f1 | sort -u > ${id_filename}
description_table_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "gene-descriptions")
cat ${id_description_table} | sort -k1 | join -a 1 -1 1 -2 1 ${id_filename} - | sed -e 's/ /|/' | tr '|' '\t' > ${description_table_filename}
regions_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "regions")
cat ${bed_filename} | cut -f4 | sort -u > ${regions_filename}
region_table_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "gene-regions")
python generate-gene-region-table.py ${bed_filename} > ${region_table_filename}

# Analysis of lost gene IDs ...
echo "Analysis of lost gene IDs ..."
lostids_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "lost-gene-ids")
cat ${ucsc_table} | sed '1d' | cut -f13 | sort -u > $outputdir/geneids-all.tmp
grep -vxF -f ${id_filename} $outputdir/geneids-all.tmp > ${lostids_filename}
rm -f $outputdir/geneids-all.tmp
nr_of_lost_genes=`cat ${lostids_filename} | wc -l`
echo "Lost ${nr_of_lost_genes}# genes ..."

# Create CisTargetX ini file ...
echo "Create CisTargetX ini file ..."
full_outputdir=`cd ${outputdir}; pwd`
full_regions_filename=`python -c "import os.path; print os.path.realpath('${regions_filename}');"`
full_region_table_filename=`python -c "import os.path; print os.path.realpath('${region_table_filename}');"`
cat cistargetx-template.ini \
		| sed -e "s!HOME!${HOME}!; s!DATADIR!${DATADIR}!; s!OUTPUTDIR!${full_outputdir}!; s!REGIONSFILE!${full_regions_filename}!; s!LUTFILE!${full_region_table_filename}!" \
		> $outputdir/cistargetx-template.ini

# Create FASTA file ...
echo "Create FASTA file ..."
fasta_filename=$(create_filename ${base_genome_id} ${outputdir} ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} "fasta")
twoBitToFa -bed=${bed_filename} ${base_genome_2bit_file} ${fasta_filename}
#Usage of -bed instead of -seqList is easier and safe ...
#awk '{printf "%s:%d-%d\n", $1, $2, $3}' ${bed_filename} | twoBitToFa ${base_genome_2bit_file} ${fasta_filename} -seqList=stdin

# Liftover procedure ...
echo "Liftover procedure ..."
${DATADIR}/genomes/generate-liftover-fasta.sh ${bed_filename} \
		${base_genome_id} $outputdir $(create_region_description ${upstream_extension_in_bp} ${delineation} ${downstream_extension_in_bp} ${intronic_extension_in_bp}) ${genomes_liftover_table}
		
# Copy ini file to output folder ...
cp ${inifile} ${outputdir}
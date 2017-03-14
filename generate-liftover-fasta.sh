#!/usr/bin/env bash

ICISTARGET_DATA_DIR='/media/data/lcb/icistarget/data';


# Parameters ...
if [ $# -ne 5 ]; then
	printf 'Wrong number of arguments.\n\nUsage: %s input_bed_file from_assembly output_dir suffix liftover_genomes_lst_file\n\n' "${0}";
	exit 2;
fi


input_bed_file="${1}";
from_assembly="${2}";
output_dir="${3}";
suffix="${4}";
liftover_genomes_lst_file="${5}";


# Check if the needed programs exist ...
if [ $(type liftOver > /dev/null 2>&1; echo $?;) -ne 0 ] ; then
    printf 'ERROR: Add "liftOver" to ${PATH}.\n';
    exit 1;
fi

if [ $(type twoBitToFa > /dev/null 2>&1; echo $?;) -ne 0 ] ; then
    printf 'ERROR: Add "twoBitToFa" to ${PATH}.\n';
    exit 1;
fi


# Performing liftOver to other species ...
to_assemblies=$(
    awk -F '\t' -v "from_assembly=${from_assembly}" '
        {
            if ($3 == from_assembly && $1 != from_assembly) {
                print $1;
            }
        }
    ' "${liftover_genomes_lst_file}";
)


for to_assembly in ${to_assemblies}; do
	# (Convert first character of to_assembly variable to uppercase with ${to_assembly^}.)
	to_assembly_first_cap="${to_assembly^}";

	chain_file="${ICISTARGET_DATA_DIR}/genomes/liftOver/${from_assembly}To${to_assembly_first_cap}.over.chain.gz";
	twobit_file="${ICISTARGET_DATA_DIR}/genomes/2bit/${to_assembly}.2bit";
	output_bed_file="${output_dir}/${to_assembly}.mapped.bed";
	log_file="${output_dir}/${to_assembly}${suffix}.liftover.log";

	printf "\nLiftover: ${from_assembly} => ${to_assembly} ...\n";
	# Relaxing stringency in parameters ...
	liftOver "${input_bed_file}" "${chain_file}" "${output_bed_file}" "${log_file}" -minMatch=0.1 -multiple || exit 1;

	printf "\nCreating FASTA file ...\n";
	twoBitToFa -bed="${output_bed_file}" "${twobit_file}" "${output_dir}/${to_assembly}${suffix}.fasta" || exit 1;
done

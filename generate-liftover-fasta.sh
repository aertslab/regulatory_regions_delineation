#!/usr/bin/env bash

# Constants ...
GENOMESDIR=$DATADIR/genomes
GENOMESLST=$GENOMESDIR/genomes.lst

#Parameters ...
if [ $# -eq 4 ]; then
	genomesfile=$GENOMESLST
elif [ $# -eq 5 ]; then
	genomesfile=$5
else
	echo "Wrong number of arguments."
	exit 2 
fi
inputfile=$1
refgenome=$2
outputdir=$3
suffix=$4

# Performing liftOver to other species ...
species=`cat $genomesfile | awk -F '\t' -v r=$refgenome '$3==r {print $1}' | grep -v $refgenome`
for curspecies in $species; do
	capital=`echo -n "${curspecies:0:1}" | tr "[:lower:]" "[:upper:]"`;
	capspecies=`echo -n "${capital}${curspecies:1}"`
	chainfile="$GENOMESDIR/${refgenome}To${capspecies}.over.chain"
	bitfile="$GENOMESDIR/$curspecies.2bit"
	bedfile="$outputdir/$curspecies.mapped.bed"
	logfile="$outputdir/$curspecies$suffix.liftover.log"
	echo "Liftover: $refgenome => $curspecies"
	# Relaxing stringency in parameters ...
	liftOver $inputfile $chainfile $bedfile $logfile -minMatch=0.1 -multiple
	echo "Creating FASTA file"
	twoBitToFa -bed=$bedfile $bitfile $outputdir/$curspecies$suffix.fasta
	#Usage of -bed instead of -seqList is easier and safe ...
	#awk '{printf "%s:%d-%d\n", $1, $2, $3}' $bedfile | twoBitToFa $bitfile $outputdir/$curspecies$suffix.fasta -seqList=stdin
	#rm -f $bedfile
done
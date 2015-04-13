#!/bin/bash
set -euo pipefail

#
#Script to generate simulated reads up to a specified coverage from
# a folder with fasta files. In particular this was used to
# simulate viral reads from the NCBI virus genomes
# 

genomes_folder="/labcommon/db/fasta/ncbi_genome/virus"
virus="Abaca_bunchy_top_virus_uid28697"

target_coverage=10
read_size=300

for virus in $(ls $genomes_folder);
do
	#Skip non folders!
	if [ ! -d ${genomes_folder}/${virus} ]; then continue; fi

	#Create dir for each virus
	mkdir -p ${virus}
	#Prefix for output files
	OUT_PREFIX=${virus}/${virus}

	#Create temp fasta file from concatenating all fna files
	cat ${genomes_folder}/${virus}/*.fna > ${OUT_PREFIX}.fa 

	#Calculate genome size
	genome_size=$(grep -vG "^>" ${OUT_PREFIX}.fa | tr -d '\n' | wc -m)
	#Calculate number of reads to reach {target_coverage} coverage (5X default)
	num_reads=$((genome_size*target_coverage/(2*read_size)))

	echo ${virus}
	echo "Genome size: ${genome_size}"
	echo "Num reads ${num_reads}"
	echo

	#Simulate reads
	#-d stands for insert size
	wgsim -1 ${read_size} -2 ${read_size} -d 500 -N ${num_reads} ${OUT_PREFIX}.fa ${OUT_PREFIX}_1.fq ${OUT_PREFIX}_2.fq

	#Delete tmp fasta file
	rm ${OUT_PREFIX}.fa
done

#!/bin/bash
#Script to get stats from libraries run with viral discovery or general metagenomics pipeline

#Requires:
# -samtools
# -assemblathon_stats.pl from https://github.com/ucdavis-bioinformatics/assemblathon2-analysis
set -euo pipefail

FOLDERS=`ls -d P*`

#Get stats from mapping to human
for x in $FOLDERS;
do
	echo $x >> human_reads_pe.txt ;
	samtools flagstat "$x"/contamination_rm/*_bwahg_pe.bam >> human_reads_pe.txt; 
	echo $x >> human_reads_se.txt ;
	samtools flagstat "$x"/contamination_rm/*_bwahg_se.bam >> human_reads_se.txt;
done

#Get assemblathon stats
ASM_STATS=asm_stats.txt
for x in $FOLDERS;
do
	echo $x >> $ASM_STATS;
	assemblathon_stats.pl "$x"/assembly/*allctgs.fa >> $ASM_STATS;
	echo "Paired ends mapped to contigs" >> $ASM_STATS
	samtools view -hSb "$x"/assembly/singletons/pe_to_contigs.sam | samtools flagstat - >> $ASM_STATS; 
	echo "Single ends mapped to contigs" >> $ASM_STATS
	samtools view -hSb "$x"/assembly/singletons/se_to_contigs.sam | samtools flagstat - >> $ASM_STATS; 
	echo >> $ASM_STATS
done

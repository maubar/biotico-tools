#!/bin/bash
set -euo pipefail

candidate_file=$1
folder=$2
out_file=$3

for asm in masurca sga fermi raymeta; 
do
	q -t "select c7 from $candidate_file where c8='$asm'" > tmp.seqids
	if [ -s tmp.seqids ];
	then
		sed -i "s/ .*//" tmp.seqids
		echo $asm
		echo `wc -l tmp.seqids`' seqs to extract'
		seqtk subseq $folder/assembly/*_asm_"$asm"_ctgs_filt.fa tmp.seqids >> $out_file
		echo `grep -G "^>" $out_file | wc -l`' seqs extracted'
		echo
	fi
done
rm tmp.seqids

#Paired-end
q -t "select c7 from $candidate_file where c8='pe'" > query.seqids
if [ -s query.seqids ];
then
	sed "s/\/[12]//" query.seqids | sort | uniq | awk '{print $0"/1";print $0"/2"}' > tmp.seqids
	echo "pe"
	echo `wc -l tmp.seqids`' seqs to extract'
	seqtk subseq $folder/tax_assign/reads_fa/*_rmcont_pe.fa tmp.seqids >> $out_file
	echo `grep -G "^>" $out_file | wc -l`' seqs extracted'
	echo
fi
rm query.seqids

#Single-end
q -t "select c7 from $candidate_file where c8='se'" > tmp.seqids
if [ -s tmp.seqids ];
then
	sed -i "s/\/[12]//" tmp.seqids
	echo "se"
	echo `wc -l tmp.seqids`' seqs to extract'
	seqtk subseq $folder/tax_assign/reads_fa/*_rmcont_se.fa tmp.seqids >> $out_file
	echo `grep -G "^>" $out_file | wc -l`' seqs extracted'
	echo
fi

rm tmp.seqids

# Viral discovery pipeline for Illumina (MiSeq) data

# Author: Mauricio Barrientos-Somarribas
# Email:  mauricio.barrientos@ki.se

# Copyright 2014 Mauricio Barrientos-Somarribas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# To-do's
#@TODO: Process single-ends from quality filtering
#@TODO: Merge overlapping paired-ends into single end- U-Search/pandaseq
#@TODO: Dataset reduction: Choose Diginorm/Fixseq/none?
#			-Try 3 pass Diginorm
#@TODO: Try other assemblers: SGA, CABOG
#@TODO: Extract singletons via mapping/file parsing after assembly
#@TODO: Classify singletons with Kraken
#@TODO: Include protein homology based classification

#Make parameters
SHELL := /bin/bash

input_files := $(wildcard reads/*.fastq.gz) $(wildcard reads/*.fq.gz) $(wildcard reads/*.fq) $(wildcard reads/*.fastq)
output_basename = P431_10x

#Run params
threads:=16

#Databases
db_path:=/labcommon/db
bwa_hg_ref:= $(db_path)/iGenomes/Homo_Sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
blastdb_folder:=$(db_path)/blastdb
kraken_db:=$(db_path)/krakendb/kraken140311/
#Mini Kraken
#kraken_db:=$(db_path)/krakendb/minikraken_20140104/

#Tool parameters
blast_params:= -evalue 10 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

#Logging info
log_file := $(output_basename)_$(shell date +%s).log
write_to_log := >( tee -a $(log_file) >&2 )

#Prefixes for output files
nesoni_pre := $(output_basename)_q20h
bwa_pre := $(nesoni_pre)_grch37
sortsam_pre := $(bwa_pre)_sort
filtersam_pre := $(sortsam_pre)_filter
sam2fq_pre := $(filtersam_pre)_sam2fq
diginorm_pre:=$(sam2fq_pre)_dgnrm

#Delete produced files if step fails
.DELETE_ON_ERROR:

all: kraken_ray.report kraken_fermi.report kraken_abyss.report

#Quality filtering
$(nesoni_pre)_R%.fq.gz: $(input_files)
	echo -e "\nNesoni Quality Filtering\n\n" > $(write_to_log)
	nesoni clip --homopolymers yes --quality 20 --length 75 \
	--out-separate yes $(nesoni_pre) pairs: $^ 2> $(write_to_log)

$(bwa_pre).sam: $(nesoni_pre)_R1.fq.gz $(nesoni_pre)_R2.fq.gz
	echo -e "\nMapping reads to GRCh37 with BWA MEM\n\n" > $(write_to_log)
	bwa mem -t $(threads) -T 20 \
	-M $(bwa_hg_ref) $^ > $@ 2> $(write_to_log)

#Convert from sam to bam removing secondary mappings
$(bwa_pre).bam: $(bwa_pre).sam
	echo -e "\nConverting .sam to .bam\n\n" > $(write_to_log)
	samtools view -F 256 -hSb -o $@ $^ 2> $(write_to_log)

#Sort files by queryname
$(sortsam_pre).bam: $(bwa_pre).bam
	echo -e "\nSort bam file by queryname\n\n" > $(write_to_log)
	run_picard SortSam.jar INPUT=$^ OUTPUT=$@ SORT_ORDER=queryname 2> $(write_to_log)

#Keep only reads that did not map confidently (with both pairs)
$(filtersam_pre).bam: $(sortsam_pre).bam
	echo -e "\nExtract reads that did not map to GRCh37\n\n" > $(write_to_log)
	run_picard FilterSamReads.jar INPUT=$^ OUTPUT=$@ FILTER=excludeAligned SORT_ORDER=queryname WRITE_READS_FILES=False 2> $(write_to_log)

#Convert unmapped reads to Fastq for assembly
#$(sam2fq_pre).R%.fq: $(filtersam_pre).bam
#	run_picard SamToFastq.jar INPUT=$^ FASTQ=$(sam2fq_pre).R1.fq SECOND_END_FASTQ=$(sam2fq_pre).R2.fq $(write_to_log)

#Alternatively use interleaved paired-end fastq format
$(sam2fq_pre).fq: $(filtersam_pre).bam
	echo -e "\nConvert bam to interleaved fastq\n\n" > $(write_to_log)
	run_picard SamToFastq.jar INPUT=$^ FASTQ=$@ INTERLEAVE=True 2> $(write_to_log)

#***********************Dataset reduction*********************************

#Apply one pass diginorm procedure, normalize to 20X cov.
#-N number of hashes 4, -x min size of hash, -k kmer size, -C coverage cutoff
$(diginorm_pre).fq: $(sam2fq_pre).fq
	echo -e "\nApply one-pass diginorm, coverage 20X\n\n" > $(write_to_log)
	normalize-by-median.py -C 20 --paired -k 20 -N 4 -x 2.5e8 $^ 2> $(write_to_log)
	mv $(sam2fq_pre).fq.keep $@

#*************************Assembly Step**********************************************
#*********************Meta-Ray with diginorm*************
raymeta/Contigs.fasta: $(diginorm_pre).fq
	echo -e "\nAssembling reads with Ray Meta\n\n" > $(write_to_log)
	mpiexec -n 16 Ray Meta -i $^ -o raymeta 2> $(write_to_log)

raymeta_contigs.fa: raymeta/Contigs.fasta
	ln -s $^ $@

#*********************Fermi with diginorm******************
#Runs Fermi assembler until the 4th step.
#Qualities are not neccessary for Kraken/Blast classification
fermi/fmdef.p4.fa.gz: $(diginorm_pre).fq
	mkdir -p fermi
	cd fermi && run-fermi.pl -t $(threads) -c ../$^ > assembly.mak 2> $(write_to_log)
	cd fermi && $(MAKE) -f assembly.mak -j $(threads) $@ 2> $(write_to_log)

fermi_contigs.fa: fermi/fmdef.p4.fa.gz
	gunzip -c $^ > $@

#*************************Abyss with diginorm***********************
abyssk47-contigs.fa: $(diginorm_pre).fq
	mkdir -p abyss
	cd abyss && abyss-pe k=47 name='abyssk47' in='$^' np=16

#*************************Run Kraken on contigs**************************
#Other flags: --fastq-input
kraken_raymeta.out kraken_fermi.out kraken_abyss.out: kraken_%.out : %_contigs.fa
	echo -e "\nClassifying $* contigs with Kraken\n\n" > $(write_to_log)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2> $(write_to_log)

kraken_raymeta.report kraken_fermi.report kraken_abyss.report: kraken_%.report : kraken_%.out
	echo -e "\nCreating Kraken report for $* \n\n" > $(write_to_log)
	kraken-report --db $(kraken_db) $^ > $@ 2> $(write_to_log)

#Run blast on contigs and output filtered results with <tool/script>
#@TODO: Create script to filter blast xml output
# megablast.%.blst: contigs.fq
# 	blastn -task megablast $(blast_params) $(megablast_params) -query $^ -db $(blastdb_path)/nt/nt > megablast.MAM.blst
#  	touch megablast.BCT.blst
# 	blastn -task megablast $(blast_params) $(megablast_params) -query $^ -db $(blastdb_path)/nt.split/nt.MAM > megablast.MAM.blst
# 	blastn -task megablast $(blast_params) $(megablast_params) -query $^ -db $(blastdb_path)/nt.split/nt.BCT > megablast.BCT.blst

#@TODO: Process blast output and output results
# results.txt: megablast.MAM.blst megablast.BCT.blst
# 	touch results.txt

.PHONY :clean
clean:
	-rm *.fq.gz *.sam *.log $(output_basename)_*
	-rm -r raymeta/ fermi/ abyss/
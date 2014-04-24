#Make parameters
SHELL := /bin/bash

input_files := $(wildcard reads/*.fastq.gz) $(wildcard reads/*.fq.gz) $(wildcard reads/*.fq) $(wildcard reads/*.fastq)
output_basename = P431_101

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

all: kraken_fermi.report kraken_ray.report kraken_fermi_nodigi.report kraken_ray_nodigi.report

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

#@TODO: Use U-Search/pandaseq to merge paired-ends with overlapping thingy



#***********************Dataset reduction - Diginorm or else *************
#@TODO: Choose between Diginorm or Fixseq?
#Apply one pass diginorm procedure, normalize to 20X cov.
#-N number of hashes 4, -x min size of hash, -k kmer size, -C coverage cutoff
$(diginorm_pre).fq: $(sam2fq_pre).fq
	echo -e "\nApply one-pass diginorm, coverage 20X\n\n" > $(write_to_log)
	normalize-by-median.py -C 20 --paired -k 20 -N 4 -x 2.5e8 $^ 2> $(write_to_log)
	mv $(sam2fq_pre).fq.keep $@
#@TODO: Try diginorm 2/3 pass
#	#Trim 3' ends with rare k-mers
# 	filter-abund.py $^ $@
#	#Renormalize coverage to 5X
# 	normalize-by-median.py -C 5 -k 20 -N 4 -x 1e8 $@


#*************************Assembly and Classification step****************

#*********************Meta-Ray with diginorm*************
ray_assembly/Contigs.fasta: $(diginorm_pre).fq
	echo -e "\nAssembling reads with Ray Meta\n\n" > $(write_to_log)
	mpiexec -n 16 Ray Meta -i $^ -o ray_assembly 2> $(write_to_log)

kraken_ray.out: ray_assembly/Contigs.fasta
	echo -e "\nClassifying Ray's scaftigs with Kraken\n\n" > $(write_to_log)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2> $(write_to_log)

kraken_ray.report: kraken_ray.out
	echo -e "\nCreating Kraken report\n\n" > $(write_to_log)
	kraken-report --db $(kraken_db) $^ > $@ 2> $(write_to_log)

#*********************Fermi with diginorm******************
#Runs Fermi assembler until the 4th step. Qualities are not neccessary for classification
fmdef.p4.fa.gz: $(diginorm_pre).fq
	run-fermi.pl -t $(threads) -c $^ > assembly.mak 2> $(write_to_log)
	$(MAKE) -f assembly.mak -j $(threads) $@ 2> $(write_to_log)

#@TODO: Run a fastq join between assembled contigs and singletons
scaftigs.fa: fmdef.p4.fa.gz
	gunzip -c $^ > $@

#Run Kraken on dataset
#Other flags: --fastq-input
kraken_fermi.out: scaftigs.fa
	echo -e "\nClassifying Fermi's scaftigs with Kraken\n\n" > $(write_to_log)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2> $(write_to_log)

kraken_fermi.report: kraken_fermi.out
	echo -e "\nCreating Kraken report\n\n" > $(write_to_log)
	kraken-report --db $(kraken_db) $^ > $@ 2> $(write_to_log)


#*********************Ray Meta - no diginorm *****************************
ray_assembly_nodigi/Contigs.fasta: $(sam2fq_pre).fq
	echo -e "\nAssembling reads with Ray - no diginorm\n\n" > $(write_to_log)
	mpiexec -n 16 Ray Meta -i $^ -o ray_assembly_nodigi 2> $(write_to_log)

kraken_ray_nodigi.out: ray_assembly_nodigi/Contigs.fasta
	echo -e "\nClassifying Ray's scaftigs with Kraken\n\n" > $(write_to_log)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2> $(write_to_log)

kraken_ray_nodigi.report: kraken_ray_nodigi.out
	echo -e "\nCreating Kraken report\n\n" > $(write_to_log)
	kraken-report --db $(kraken_db) $^ > $@ 2> $(write_to_log)

#*********************Fermi no diginorm*****************
fermi_nodigi.p4.fa.gz: $(sam2fq_pre).fq
	run-fermi.pl -p fermi_nodigi -t $(threads) -c $^ 2> assembly_nodigi.mak $(write_to_log)
	$(MAKE) -f assembly_nodigi.mak -j $(threads) $@ 2> $(write_to_log)

scaftigs_nodigi.fa: fermi_nodigi.p4.fa.gz
	gunzip -c $^ > $@ 2> $(write_to_log)

kraken_fermi_nodigi.out: scaftigs_nodigi.fa
	echo -e "\nClassifying Fermi's scaftigs with Kraken\n\n" > $(write_to_log)
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ 2> $(write_to_log)

kraken_fermi_nodigi.report: kraken_fermi_nodigi.out
	echo -e "\nCreating Kraken report\n\n" > $(write_to_log)
	kraken-report --db $(kraken_db) $^ > $@ 2> $(write_to_log)


#@TODO: SGA Assembly

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
	-rm -r ray_assembly

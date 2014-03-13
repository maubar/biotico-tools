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
minikraken_db:=$(db_path)/krakendb/minikraken_20140104/

#Tool parameters
blast_params:= -evalue 10 -num_threads $(threads) -max_target_seqs 10 -outfmt 5 -show_gis
megablast_params:= -reward 2 -penalty -3 -gapopen 5 -gapextend 2
blastn_params:= -reward 4 -penalty -5 -gapopen 12 -gapextend 8

#Logging info
log_file := $(output_basename)_$(shell date +%s).log
write_to_log := 2> >( tee -a $(log_file) >&2 )

#Prefixes for output files
nesoni_pre := $(output_basename)_q20h
bwa_pre := $(nesoni_pre)_grch37
sortsam_pre := $(bwa_pre)_sort
filtersam_pre := $(sortsam_pre)_filter
sam2fq_pre := $(filtersam_pre)_sam2fq
diginorm_pre:=$(sam2fq_pre)_dgnrm

#Delete produced files if step fails
.DELETE_ON_ERROR:

all: kraken.report

#Quality filtering
$(nesoni_pre)_R%.fq.gz: $(input_files)
	nesoni clip --homopolymers yes --quality 20 --length 75 \
	--out-separate yes $(nesoni_pre) pairs: $^ $(write_to_log)

$(bwa_pre).sam: $(nesoni_pre)_R1.fq.gz $(nesoni_pre)_R2.fq.gz
	bwa mem -t $(threads) -T 20 \
	-M $(bwa_hg_ref) $^ > $@ $(write_to_log)

#Convert from sam to bam removing secondary mappings
$(bwa_pre).bam: $(bwa_pre).sam
	samtools view -F 256 -hSb -o $@ $^ $(write_to_log)

#Sort files by queryname
$(sortsam_pre).bam: $(bwa_pre).bam
	run_picard SortSam.jar INPUT=$^ OUTPUT=$@ SORT_ORDER=queryname $(write_to_log)

#Keep only reads that did not map confidently (with both pairs)
$(filtersam_pre).bam: $(sortsam_pre).bam
	run_picard FilterSamReads.jar INPUT=$^ OUTPUT=$@ FILTER=excludeAligned SORT_ORDER=queryname WRITE_READS_FILES=False $(write_to_log)

#Convert unmapped reads to Fastq for assembly
#$(sam2fq_pre).R%.fq: $(filtersam_pre).bam
#	run_picard SamToFastq.jar INPUT=$^ FASTQ=$(sam2fq_pre).R1.fq SECOND_END_FASTQ=$(sam2fq_pre).R2.fq $(write_to_log)

#Alternatively use interleaved paired-end fastq format
$(sam2fq_pre).fq: $(filtersam_pre).bam
	run_picard SamToFastq.jar INPUT=$^ FASTQ=$@ INTERLEAVE=True $(write_to_log)

#@TODO: Use U-Search to merge paired-ends with overlapping thingy


#@TODO: Choose between Diginorm or Fixseq?
#Apply one pass diginorm procedure, normalize to 20X cov.
#-N number of hashes 4, -x min size of hash, -k kmer size, -C coverage cutoff
$(diginorm_pre).fq: $(sam2fq_pre).fq
	normalize-by-median.py -C 20 --paired -k 20 -N 4 -x 2.5e8 $^ $(write_to_log)
	mv $(sam2fq_pre).fq.keep $(diginorm_pre).fq
#@TODO: Try diginorm 2/3 pass
#	#Trim 3' ends with rare k-mers 
# 	filter-abund.py $^ $@
#	#Renormalize coverage to 5X
# 	normalize-by-median.py -C 5 -k 20 -N 4 -x 1e8 $@


#@TODO: Try other assemblers (Velvet, MetaRay, MetaIDBA(now IDBA-UC :S) )

#Assembly with MetaRay
#@TODO: Choose Scaffolds or Contigs?
ray_assembly/Scaffolds.fasta: $(diginorm_pre).fq
	mpiexec -n 16 Ray Meta -i $^ -o ray_assembly $(write_to_log)

#Assembly with Fermi
#fmdef.p5.fq.gz: $(sam2fq_pre).fq.keep
#	#run-fermi.pl -t $(threads) -c $^ > assembly.mak $(write_to_log)
#	run-fermi.pl -t 8 -c $^ > assembly.mak $(write_to_log)
#	$(MAKE) -f assembly.mak -j $(threads) $(write_to_log)

#@TODO: Run a fastq join between assembled contigs and singletons
#contigs.fq: fmdef.p5.fq.gz
#	gunzip -c $^ > contigs.fq

#Run Kraken on dataset
#Other flags: --fastq-input 
kraken.out: ray_assembly/Scaffolds.fasta
	kraken --db $(kraken_db) --threads $(threads) $^ > $@ $(write_to_log)

kraken.report: kraken.out
	kraken-report --db $(kraken_db) $^ > $@ $(write_to_log)

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

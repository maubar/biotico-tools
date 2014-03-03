#Make parameters
SHELL := /bin/bash

input_files := $(wildcard reads/*.fastq.gz)
output_basename = test

#Databases
bwa_hg_ref := /proj/b2011088/db/GRCh37/iGenome/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
blastdb_folder:= /proj/b2011088/nobackup

#Tool parameters
blast_params:= -evalue 10 -num_threads 16 -max_target_seqs 10 -outfmt 5 -show_gis
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

all: results.txt

#Quality filtering
$(nesoni_pre)_R%.fq.gz: $(input_files)
	nesoni clip --homopolymers yes --quality 20 --length 75 \
	--out-separate yes $(nesoni_pre) pairs: $^ $(write_to_log)

$(bwa_pre).sam: $(nesoni_pre)_R1.fq.gz $(nesoni_pre)_R2.fq.gz
	bwa mem -t 8 -T 20 \
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
$(sam2fq_pre).R%.fq: $(filtersam_pre).bam
	run_picard SamToFastq.jar INPUT=$^ FASTQ=$(sam2fq_pre).R1.fq SECOND_END_FASTQ=$(sam2fq_pre).R2.fq $(write_to_log)

#@TODO: Diginorm step before assembly (debrujin graph)

#@TODO: Try other assemblers (Velvet, MetaRay, MetaIDBA)
#Run assembly with Fermi
fmdef.p5.fq.gz: $(sam2fq_pre).R1.fq $(sam2fq_pre).R2.fq 
	run-fermi.pl -t 16 -P $^ > assembly.mak $(write_to_log)
	$(MAKE) -f assembly.mak -j 16 $(write_to_log)

#@TODO: Run a fastq join between assembled contigs and singletons
contigs.fq: fmdef.p5.fq.gz
	gunzip -c $^ > contigs.fq

#Run blast on contigs and output filtered results with <tool/script>
#@TODO: Create script to filter blast xml output
megablast.%.blst: contigs.fq
	blastn -task megablast $(blast_params) $(megablast_params) -query $^ -db $(blastdb_path)/nt.split/nt.MAM > megablast.MAM.blst
	blastn -task megablast $(blast_params) $(megablast_params) -query $^ -db $(blastdb_path)/nt.split/nt.BCT > megablast.BCT.blst

#@TODO: Process blast output and output results
results.txt: megablast.MAM.blst megablast.BCT.blst
	touch results.txt

.PHONY :clean-all
clean-all:
	-rm *.fq.gz *.sam

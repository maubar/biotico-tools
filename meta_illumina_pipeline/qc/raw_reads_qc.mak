#Make parameters
SHELL := /bin/bash

R1:= $(wildcard *R1*.fastq.gz *_1.fastq.gz)
R2:= $(wildcard *R2*.fastq.gz *_2.fastq.gz)

#Logging
log_name := qc_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

#Run params
threads:=16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: sga_preproc.preqc k17.hist.pdf $(basename $(R1))_fastqc.zip $(basename $(R2))_fastqc.zip

#*************************************************************************
#SGA PREQC
#*************************************************************************
#PREQC: First, preprocess the data to remove ambiguous basecalls
sga_preproc.fq: $(R1) $(R2)
	sga preprocess --pe-mode 1 -o $@ $^ >&2 2>> $(log_file)

# PREQC: Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
sga_preproc.sai: sga_preproc.fq
	sga index -a ropebwt -t $(threads) --no-reverse $(notdir $^) >&2 2>> $(log_file)

#Run SGA preqc
%.preqc: %.fq %.sai
	sga preqc -t $(threads) --force-EM $< > $@ 2>> $(log_file)
	cd $(dir $@) && sga-preqc-report.py $(notdir $@)

#*************************************************************************
#FASTQC
#*************************************************************************
%.fastq_fastqc.zip: %.fastq.gz
	fastqc --noextract -k 10 $^ 2>> $(log_file)

#*************************************************************************
#JELLYFISH 2
#*************************************************************************
k17.jf: sga_preproc.fq
	jellyfish2 count -s 8G -C -m 17 -t $(threads) -o $@ $^ 2>> $(log_file)

k17.hist: k17.jf
	jellyfish2 histo -t $(threads) $^ -o $@ 2>> $(log_file)

k17.hist.pdf: k17.hist
	Rscript plot_kmer_histogram.R $^ $@

#*************************************************************************
#CLEANING RULES
#*************************************************************************
.PHONY: clean-tmp clean-out

clean-tmp:
	-rm sga_preproc.{fq,sai}
	-rm *.log #Makefile log
	-rm *.sai *.bwt
	-rm *.jf

clean-out:
	-rm *.fq_fastqc.zip #Fastqc
	-rm *.hist *.hist.pdf #Jellyfish
	-rm *.preqc *.pdf #SGA preqc

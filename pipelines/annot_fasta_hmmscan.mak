SHELL:=/bin/bash

in_fasta := $(wildcard input/*.fa)

SAMPLE:= sample
threads := 16

PfamA_hmm := /labcommon/db/hmmerdb/Pfam-A/Pfam-A.hmm
vFamA_hmm := /labcommon/db/hmmerdb/vFam-A/vFam-A_2014.hmm

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: hmmscan/$(SAMPLE)_emboss_orfs_hmmscan_PfamA.out
all: hmmscan/$(SAMPLE)_emboss_orfs_hmmscan_vFamA.out

#*********************************************************
# Predict ORFs with EMBOSS getorf
#*********************************************************
orf/$(SAMPLE)_emboss_find0.fa: $(in_fasta)
	mkdir -p $(dir $@)
	getorf -sequence $< -outseq $@ -find 0

orf/$(SAMPLE)_emboss_find1.fa: $(in_fasta)
	mkdir -p $(dir $@)
	getorf -sequence $< -outseq $@ -find 1

orf/$(SAMPLE)_emboss_orfs.fa: orf/$(SAMPLE)_emboss_find0.fa orf/$(SAMPLE)_emboss_find1.fa
	cat  <( awk '/^>/ {print gensub(/(_[0-9]+) (\[.+\])/,"\\1_find0 \\2",$$0); } ! /^>/{print $$0;}' $< ) \
	 	 <( awk '/^>/ {print gensub(/(_[0-9]+) (\[.+\])/,"\\1_find1 \\2",$$0); } ! /^>/{print $$0;}' $(word 2,$^) ) > $@

#*********************************************************
# Run hmmscan Pfam-A against the orfs
#*********************************************************
hmmscan_opts = --noali --tblout $(basename $@).tbl --domtblout $(basename $@).domtbl

hmmscan/%_hmmscan_PfamA.out: orf/%.fa | $(PfamA_hmm)
	mkdir -p $(dir $@)
	hmmscan --cpu $(threads) $(hmmscan_opts) -o $@ $| $<

hmmscan/%_hmmscan_vFamA.out: orf/%.fa | $(vFamA_hmm)
	mkdir -p $(dir $@)
	hmmscan --cpu $(threads) $(hmmscan_opts) -o $@ $| $<

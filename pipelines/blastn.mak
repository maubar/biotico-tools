SHELL:=/bin/bash

ifndef out_prefix
$(error Variable 'out_prefix' is not defined)
endif

ifndef in_fasta
$(error Variable 'in_fasta' is not defined)
endif

ifndef blastdb
$(error Variable 'blastdb' is not defined)
endif

threads := 16

.DELETE_ON_ERROR:

.SECONDARY:

.PHONY: all

all: $(out_prefix)/$(out_prefix).tsv

#Indexes

#****************************************************
#Blast parameters
#***************************************************

blast_params:=-task blastn -num_threads $(threads) -max_target_seqs 10 -outfmt 11

$(out_prefix)/$(out_prefix).asn: $(in_fasta)
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb) -query $< -out $@

#*******************************************************
# Formatters
#*******************************************************
%.html: %.asn
	blast_formatter -archive $< -html -out $@

blast_tsv_cols := qseqid sseqid pident qcovs score qlen qstart qend slen sstart send length mismatch gapopen

%.tsv: %.asn
	blast_formatter -archive $< -outfmt '6 $(blast_tsv_cols)' -out $@
	echo "$(blast_tsv_cols)" > $(dir $@)/blast_tsv_columns.txt

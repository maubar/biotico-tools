#Make parameters
SHELL := /bin/bash

#SGA parameters
sga_ec_kmer := 41
sga_cov_filter := 2


#*************************************************************************
#SGA quality filtering steps - not used
#*************************************************************************
#1) Preprocess the data to remove ambiguous basecalls
sga_qc/$(sample_name)_sga.fq: $(R1) $(R2)
	$(SGA_BIN) preprocess --pe-mode 1 -o $@ $^

#2) Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
%_sga.sai: %_sga.fq
	cd $(dir $@) && $(SGA_BIN) index -a ropebwt -t $(threads) --no-reverse $(notdir $^)

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.fq: %.k$(sga_ec_kmer).ec.fq : %.fq %.sai
	$(SGA_BIN) correct -k $(sga_ec_kmer) --discard --learn -t $(threads) -o $@ $<

# Index the corrected data
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.sai: %.k$(sga_ec_kmer).ec.sai : %.k$(sga_ec_kmer).ec.fq
	$(SGA_BIN) index -a ropebwt -t $(threads) $^

# Remove exact-match duplicates and reads with low-frequency k-mers
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.filter.pass.fa: $(sga_1).k$(sga_ec_kmer).ec.fq $(sga_1).k$(sga_ec_kmer).ec.sai
	$(SGA_BIN) filter -x $(sga_cov_filter) -t $(threads) --homopolymer-check \
		--low-complexity-check $<

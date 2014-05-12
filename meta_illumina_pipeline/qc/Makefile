#Make parameters
SHELL := /bin/bash

R1:= $(wildcard ../reads/*R1*.fastq.gz)
R2:= $(wildcard ../reads/*R2*.fastq.gz)

output_basename := P431_102_QC

#Logging
log_file := $(output_basename)_$(shell date +%s).log
write_to_log := >( tee -a $(log_file) >&2 )

#Run params
threads:=16
adaptor_file:= adapters/truseq_custom.fa
TRIMMOMATIC:= java -jar /labcommon/tools/Trimmomatic-0.32/trimmomatic-0.32.jar
sga_ec_kmer := 41
sga_cov_filter := 2

#Output basenames
nesoni_1 := nesoni_def
nesoni_2 := nesoni_custom
nesoni_3 := nesoni_match15mm2
nesoni_4 := nesoni_match15mm2_truseqht

trimmy_1 := trimmy
trimmy_nesoni := nesoni_noclip_trimmy

sga_1:= sga

nesoni_out_prefix = $(dir $@)P431_102_qf
trimmomatic_out = $(addprefix, $(dir $@)P431_102_qf_ , $(addsuffix .fq.gz R1 R1_single R2 R2_single ))

prinseq := prinseq_1

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

# all: $(nesoni_1)_R1.fq_fastqc.zip $(nesoni_2)_R1.fq_fastqc.zip $(nesoni_3)_R1.fq_fastqc.zip $(nesoni_4)_R1.fq_fastqc.zip $(trimmy_1)_R1.fq_fastqc.zip
all: $(trimmy_nesoni)_R1.fq_fastqc.zip $(trimmy_nesoni)_R2.fq_fastqc.zip


$(nesoni_1) $(nesoni_2) $(nesoni_3) $(nesoni_4) $(trimmy_1) $(trimmy_nesoni) $(prinseq):
	mkdir -p $@

$(nesoni_1)/R%.fq.gz: $(R1) $(R2) | $(nesoni_1)
	nesoni clip --length 75 --out-separate yes $(nesoni_out_prefix) pairs: $^ 2> $(write_to_log)

$(nesoni_2)/R%.fq.gz: $(R1) $(R2) | $(nesoni_2)
	mkdir $(dir $@)
	nesoni clip --adaptor-file $(adaptor_file) --length 75 --out-separate yes $(nesoni_out_prefix) pairs: $^ 2> $(write_to_log)

$(nesoni_3)/R%.fq.gz: $(R1) $(R2) | $(nesoni_3)
	nesoni clip --match 15 --max-errors 2 --length 75 --out-separate yes $(nesoni_out_prefix) pairs: $^ 2> $(write_to_log)

$(nesoni_4)/R%.fq.gz: $(R1) $(R2) | $(nesoni_4)
	nesoni clip --adaptor-file $(adaptor_file) --match 15 --max-errors 2 --length 75 --out-separate yes $(nesoni_out_prefix) pairs: $^ 2> $(write_to_log)

$(trimmy_1)/R%.fq.gz: $(R1) $(R2) | $(trimmy_1)
	$(TRIMMOMATIC) PE -threads 16 -phred33 -trimlog trimmomatic.log $^ $(trimmomatic_out) ILLUMINACLIP:$(adaptor_file):2:30:10 MINLEN:75 2> $(write_to_log)

$(trimmy_nesoni)/R%.fq.gz: $(R1) $(R2) | $(trimmy_nesoni)
	nesoni clip --adaptor-clip no --homopolymers yes --quality 20 --length 75 --out-separate yes $(dir $@)nesoni_temp pairs: $^ 2> $(write_to_log)
	$(TRIMMOMATIC) PE -threads 16 -phred33 -trimlog trimmomatic.log nesoni_temp_R1.fq.gz nesoni_temp_R2.fq.gz $(trimmy_nesoni)_R1.fq.gz $(trimmy_nesoni)_singleR1.fq.gz $(trimmy_nesoni)_R2.fq.gz $(trimmy_nesoni)_singleR2.fq.gz ILLUMINACLIP:$(adaptor_file):2:30:10 MINLEN:75 2> $(write_to_log)
	-rm nesoni_temp*

$(prinseq)/cosillo.fq: $(R1) $(R2) | $(prinseq)
	prinseq-lite.pl -fastq $(R1) -fastq2 $(R2) -params prinseq_params.txt

#Run SGA steps for qc
# First, preprocess the data to remove ambiguous basecalls
%_sga.fq: %_R1_001.fastq.gz %_R2_001.fastq.gz
	sga preprocess --pe-mode 1 -o $@ $^

# Build the index that will be used for error correction
# Error corrector does not require the reverse BWT
sga_qc/$(sga_1)_sga.sai: %_sga.sai : %_sga.fq
	sga index -a ropebwt -t $(threads) --no-reverse $^

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.fq: %.k($sga_ec_kmer).ec.fq : %.fq %.sai
	sga correct -k $(sga_ec_kmer) --discard --learn -t $(threads) -o $@ $< 

# Index the corrected data
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.sai: %.k$(sga_ec_kmer).ec.sai : %.k$(sga_ec_kmer).ec.fq
	sga index -a ropebwt -t $(threads) $^

# Remove exact-match duplicates and reads with low-frequency k-mers
sga_qc/$(sga_1).k$(sga_ec_kmer).ec.filter.pass.fa: $(sga_1).k$(sga_ec_kmer).ec.fq $(sga_1).k$(sga_ec_kmer).ec.sai
	sga filter -x $(sga_cov_filter) -t $(threads) --homopolymer-check --low-complexity-check $<

#############QUALITY VERIFICATION########################

###########FASTQC############
%.fq_fastqc.zip: %.fq.gz
	fastqc --noextract -k 10 $^

%.fq_fastqc.zip: %.fq
	fastqc --noextract -k 10 $^

##########SGA PREQC###########
#Depends on the rules for creating the index
#Run SGA preqc
%.preqc: %.fq %.sai
	sga preqc -t 16 $< > $@
	sga-preqc-report.py $@

##############JELLYFISH############

%_jf_17mer/k17_000: %.fq
	mkdir -p $*_jf_17mer
	jellyfish count -s 8G -C -m 17 -t 16 -o $(dir $@)k17 $^

%_jf_17mer/k17.jf: %_jf_17mer/k17_000
	jellyfish merge -o $@ $(dir $<)k17_*

%.k17.hist: %_jf_17mer/k17.jf
	jellyfish histo -t 16 -f $^ > $@

%.k31.hist: %.fq
	mkdir -p $*_jf_31mer
	jellyfish count -s 8G -C -m 31 -t 16 -o $*_jf_31mer/k31 $^
	cd $*_jf_31mer && jellyfish merge -o k31.jf k31_*
	jellyfish histo -t 16 -f $*_jf_31mer/k31.jf > $*.k31.hist

.PHONY: clean
clean:
	-rm *.fq.gz 
	-rm *.fq_fastqc.zip #Fastqc
	-rm *.log #Makefile log
	-rm *.log.txt #Nesoni log

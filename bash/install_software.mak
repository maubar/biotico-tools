#Requires:
# wget
# git
# GNU sed as sed
# Python 2.7 and 3, with pip installed for both
# tar , unzip
# make and cmake(for bamtools)

SHELL := /bin/bash

ROOT_FOLDER := $(shell pwd)

LINK_TO_BIN = ln -st bin/ $(1)
GET_EXECUTABLES_FROM_FOLDER=$(shell for x in $$(find $(1) -maxdepth 1 -executable -not -type d); do basename $$x;done)
MAKE_NOFLAGS := make MAKEFLAGS= 

#Python configuration.
# Commands like pip or python might require sudo or a virtualenv
pip_python2 := pip
pip_python3 := source activate py3k && pip
python2_bin := python
python3_bin := source activate py3k && python3

.PHONY: maars viral-discovery
.PHONY: fastqc jellyfish2
.PHONY: cutadapt nesoni
.PHONY: raymeta masurca

#******************************************************
#		MAARS pipeline tools
#*****************************************************
maars: seqtk samtools prinseq
maars: fastqc
maars: cutadapt nesoni
maars: bwa picard-tools
maars: megahit
maars: diamond kraken

#******************************************************
#		Viral discovery tools
#*****************************************************
viral-discovery: seqtk samtools prinseq sga
viral-discovery: fastqc
viral-discovery: cutadapt nesoni
viral-discovery: bwa picard-tools
viral-discovery: spades iva fermi raymeta masurca
viral-discovery: ncbi-blast hmmer

#Create bin folder
$(shell mkdir -p bin/)

#**************************************************************************************
#******************        General purpose tools      *********************************
#**************************************************************************************
seqtk/:
	git clone https://github.com/lh3/seqtk.git
	cd seqtk && $(MAKE_NOFLAGS)
	$(call LINK_TO_BIN,`pwd`/seqtk/seqtk)

tabtk:
	git clone https://github.com/lh3/tabtk.git
	cd tabtk && $(MAKE_NOFLAGS)
	$(call LINK_TO_BIN,`pwd`/tabtk/tabtk)

samtools/1.2: VERSION=1.2
samtools/1.2:
	wget -N https://github.com/samtools/samtools/releases/download/$(VERSION)/samtools-$(VERSION).tar.bz2
	tar -xjf samtools-*.tar.bz2
	mkdir -p samtools && mv samtools-*/ $@
	cd $@ && make

picard-tools/1.135: VERSION=1.135
picard-tools/1.135:
	wget -N https://github.com/broadinstitute/picard/releases/download/$(VERSION)/picard-tools-$(VERSION).zip
	unzip picard-tools-*.zip
	mkdir picard-tools && mv picard-tools-*/ $@
	rm picard-tools-*.zip

.PHONY: bedtools2
bedtools2:
	#wget -N https://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.23.0.tar.gz
	#tar -xzf bedtools-*.tar.gz
	cd $@ && $(MAKE_NOFLAGS)
	$(call LINK_TO_BIN,`pwd`/bedtools2/*)
	
#**************************************************************************************
#******************        QC and QF       ********************************************
#**************************************************************************************
FastQC/0.11.3: VERSION=0.11.3
FastQC/0.11.3:
	-wget -N http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v$(VERSION).zip
	mkdir -p FastQC
	unzip -d FastQC fastqc_*.zip
	mv FastQC/FastQC $@
	chmod ug+x $@/fastqc
	@echo "This installation assumes you have at least 16 gigs of ram"
	sed -i "s/Xmx250m/Xmx16g/" $@/fastqc
	rm fastqc_v*.zip

prinseq/:
	wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download -O prinseq.tar.gz
	tar -xzf prinseq.tar.gz
	mv prinseq-* prinseq/
	chmod ug+x prinseq/prinseq-lite.pl
	$(call LINK_TO_BIN,`pwd`/prinseq/prinseq-lite.pl)

jellyfish2/2.2.3: VERSION=2.2.3
jellyfish2/2.2.3: 
	-wget -N https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-$(VERSION).tar.gz
	tar -xzf jellyfish-$(VERSION).tar.gz
	mkdir -p jellyfish2 
	cd jellyfish-2*/ && ./configure --prefix=$(shell pwd)/$@ && make && make install
	rm -r jellyfish-$(VERSION) jellyfish-$(VERSION).tar.gz

cutadapt:
	$(pip_python2) install cutadapt

nesoni:
	$(pip_python2) install nesoni

fqtrim:
	-wget -N http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.94.tar.gz
	tar -xzf fqtrim-*.tar.gz
	mv fqtrim-*/ fqtrim
	cd fqtrim && make release
	$(call LINK_TO_BIN,`pwd`/fqtrim/fqtrim)

pandaseq:
	echo "Requires zlib, bzip2 and libtool"
	#sudo apt-get install zlib1g-dev libbz2-dev libltdl-dev libtool
	git clone http://github.com/neufeld/pandaseq.git/
	cd pandaseq && bash autogen.sh && ./configure --prefix=`pwd`/dist && make && make install && make clean
	$(call LINK_TO_BIN,`pwd`/pandaseq/dist/bin/*)

FLASH:
	wget -N http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz
	tar -xzf FLASH-*.tar.gz
	mv FLASH-*/ FLASH/
	cd FLASH && make
	$(call LINK_TO_BIN,`pwd`/FLASH/flash)


#**************************************************************************************
#******************        READ MAPPERS/ALIGNERS     **********************************
#**************************************************************************************
bwa/0.7.12: VERSION=0.7.12
bwa/0.7.12:
	wget -N http://downloads.sourceforge.net/project/bio-bwa/bwa-$(VERSION).tar.bz2
	tar -xjf bwa-*.tar.bz2
	mkdir -p bwa
	mv bwa-*/ $@
	cd $@ && make
	rm bwa-*.tar.bz2

bowtie2/2.2.5: VERSION=2.2.5
bowtie2/2.2.5:
	mkdir -p bowtie2
	-wget -N bowtie http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$(VERSION)/bowtie2-2.2.5-linux-x86_64.zip
	unzip bowtie2-*.zip
	mv bowtie2-*/ $@
	rm bowtie2-*.zip

#Requires python2.x and python2.x-config, where x = 6 or 7
#This is breaking, maybe requires python-dev package?
stampy/:
	wget -N http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
	tar -xzf Stampy-latest.tgz
	mv stampy-*/ stampy
	cd stampy && make
	$(call LINK_TO_BIN,`pwd`/stampy/stampy.py)

vsearch:
	wget -N https://github.com/torognes/vsearch/releases/download/v1.1.3/vsearch-1.1.3-linux-x86_64
	mv vsearch-* bin/vsearch
	chmod ug+x bin/vsearch

last:
	wget -N http://last.cbrc.jp/last-572.zip
	unzip last-572 
	mv last-*/ last
	cd $@ && make
	cd $@ && make install prefix=../bin
	 
#**************************************************************************************
#******************        ASSEMBLERS     *********************************************
#**************************************************************************************

fermi:
	wget -N https://github.com/downloads/lh3/fermi/fermi-1.1.tar.bz2
	tar -xjf fermi-*
	mv fermi-*/ fermi
	cd fermi && make
	$(call LINK_TO_BIN,`pwd`/$@/{fermi,run_fermi.pl})

megahit:
	git clone https://github.com/voutcn/megahit.git
	cd megahit && make
	$(call LINK_TO_BIN,`pwd`/$@/{megahit*,sdbg_builder_cpu})

SPAdes:
	wget -N http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
	tar -xzf SPAdes*.tar.gz && mv SPAdes-*/ SPAdes
	$(call LINK_TO_BIN,`pwd`/$@/bin/*)

#Assumes python3 installed from anaconda in a virtualenv called py3k
Fastaq:
	git clone https://github.com/sanger-pathogens/Fastaq.git
	cd Fastaq && $(python3_bin) setup.py install

kmc:
	mkdir -p kmc/
	cd kmc && wget -N http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc
	cd kmc && wget -N http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc_dump
	$(call LINK_TO_BIN,`pwd`/$@/*)

#Assumes python 3 install from anaconda in a virtual env called py3k
iva: Fastaq MUMmer smalt kmc
	$(pip_python3) install networkx
	$(pip_python3) install pysam
	wget -N https://github.com/sanger-pathogens/iva/archive/v0.11.0.tar.gz
	tar -xzf v0.11.0.tar.gz
	mv iva*/ iva
	cd iva && $(python3_bin) setup.py install

MUMmer:
	wget -N http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
	tar -xzf MUMmer*.tar.gz
	mv MUMmer*/ MUMmer
	cd MUMmer && make install
	find MUMmer/ -maxdepth 1 -executable -not -type d -exec ln -st bin/ `pwd`/{} \;

smalt:
	wget -N http://downloads.sourceforge.net/project/smalt/smalt-0.7.6-static.tar.gz
	tar -xzf smalt*.tar.gz
	mv smalt*/ smalt/
	cd smalt && mkdir -p dist && ./configure --prefix=`pwd`/dist && make install
	$(call LINK_TO_BIN,`pwd`/$@/dist/bin/*)

#SGA and sga-dependencies
sga: sparsehash bamtools jemalloc
	git clone https://github.com/jts/sga.git
	cd sga/src && ./autogen.sh
	cd sga/src && ./configure --with-sparsehash=$(ROOT_FOLDER)/sparsehash --with-bamtools=$(ROOT_FOLDER)/bamtools --with-jemalloc=$(ROOT_FOLDER)/jemalloc/lib --prefix=$(ROOT_FOLDER)
	cd sga/src && make && make install
	#Copy preqc report script to bin
	$(call LINK_TO_BIN,`pwd`/$@/src/bin/*)

sparsehash:
	wget -N https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz
	tar -xzf sparsehash*.tar.gz
	cd sparsehash-*/ && ./configure --prefix=$(ROOT_FOLDER)/sparsehash && make && make install

bamtools:
	git clone https://github.com/pezmaster31/bamtools.git
	mkdir -p bamtools/build
	cd bamtools/build && cmake ..
	cd bamtools/build && make

jemalloc:
	wget -N http://www.canonware.com/download/jemalloc/jemalloc-3.6.0.tar.bz2
	tar -xjf jemalloc-*.tar.bz2
	cd jemalloc-*/ && ./configure --prefix=$(ROOT_FOLDER)/jemalloc && make && make install

#**************************************************************************************
#******************        ASSEMBLY EVAL     *********************************************
#**************************************************************************************
assemblathon2-analysis:
	git clone https://github.com/ucdavis-bioinformatics/assemblathon2-analysis.git
	chmod ug+x $@/assemblathon-stats.pl
	$(call LINK_TO_BIN,`pwd`/$@/assemblathon-stats.pl)

#************************************************************************
#******************    GENE PREDICTION       ****************************
#************************************************************************
FragGeneScan:
	wget http://sourceforge.net/projects/fraggenescan/files/latest/download -O fgs.tar.gz
	tar -xzf fgs.tar.gz
	mv FragGeneScan*/ FragGeneScan
	cd FragGeneScan && make clean && make fgs

MetaGeneMark:
	@echo "This programs has to be manually downloaded"
	wget -N http://topaz.gatech.edu/GeneMark/tmp/GMtool_WqFda/MetaGeneMark_linux_64.tar.gz
	wget -N http://topaz.gatech.edu/GeneMark/tmp/GMtool_WqFda/gm_key_64.gz

#************************************************************************
#******************        ANNOTATION       ****************************
#************************************************************************
prokka/1.11: VERSION=1.11
prokka/1.11:
	@echo "Prokka requires the installation of Perl's libdatetime and libxml-simple"
	-wget -N http://www.vicbioinformatics.com/prokka-$(VERSION).tar.gz
	tar -xzf prokka-$(VERSION).tar.gz
	mkdir -p prokka && mv prokka-*/ $@
	$@/bin/prokka --setupdb
	rm -r prokka-$(VERSION).tar.gz

#************************************************************************
#******************        DB SEARCHES       ****************************
#************************************************************************
ncbi-blast/2.2.31+: VERSION=2.2.31
ncbi-blast/2.2.31+:
	wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(VERSION)/ncbi-blast-*-x64-linux.tar.gz
	tar -xzf ncbi-blast*.tar.gz
	mkdir -p ncbi-blast && mv ncbi-blast-*/ $@
	rm ncbi-blast-*.tar.gz

hmmer:
	wget -N ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-*-linux-intel-x86_64.tar.gz
	tar -xzf hmmer-*-linux-intel-*.tar.gz
	mv hmmer-*/ hmmer
	$(call LINK_TO_BIN,`pwd`/$@/binaries/?hmm*)
	$(call LINK_TO_BIN,`pwd`/$@/binaries/hmm*)

diamond/0.7.9: VERSION=0.7.9
diamond/0.7.9:
	wget -N http://github.com/bbuchfink/diamond/releases/download/v$(VERSION)/diamond-linux64.tar.gz
	mkdir -p $@ && cd $@ && tar -xzf ../../diamond*.tar.gz
	rm diamond*.tar.gz

kraken/0.10.5-beta: VERSION=0.10.5-beta
kraken/0.10.5-beta:
	wget -N http://ccb.jhu.edu/software/kraken/dl/kraken-$(VERSION).tgz
	tar -xzf kraken-$(VERSION).tgz
	mkdir -p $@
	cd kraken-$(VERSION)/ && bash install_kraken.sh ../$@
	rm -rf kraken-$(VERSION)/ kraken-*.tgz

#************************************************************************
#****************      VARIANT CALLING       ****************************
#************************************************************************
freebayes:
	git clone --recursive git://github.com/ekg/freebayes.git
	cd freebayes && $(MAKE_NOFLAGS)
	$(call LINK_TO_BIN,`pwd`/$@/bin/*)

vcflib:
	@echo "If rule fails, run without -r flag"
	git clone --recursive https://github.com/ekg/vcflib.git
	cd vcflib && env && $(MAKE_NOFLAGS) 
	$(call LINK_TO_BIN,`pwd`/$@/bin/*)

#************************************************************************
#****************      READ SIMULATION       ****************************
#************************************************************************
.PHONY: BEAR 
BEAR: DRISEE
	#git clone https://github.com/sej917/BEAR.git
	$(call LINK_TO_BIN,`pwd`/$@/scripts/*)

#Requires Qiime, Biopython , perl 
DRISEE:
	git clone https://github.com/MG-RAST/DRISEE.git
	$(call LINK_TO_BIN,$(addprefix $(ROOT_DIR)/$@/,$(call GET_EXECUTABLES_FROM_FOLDER,$@)))

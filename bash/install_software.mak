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
viral-discovery: seqtk samtools 
viral-discovery: fastqc
viral-discovery: cutadapt nesoni
viral-discovery: bwa picard-tools
viral-discovery: spades megahit fermi iva
viral-discovery: ncbi-blast hmmer

#Create bin folder
$(shell mkdir -p bin/)

#**************************************************************************************
#******************        General purpose tools      *********************************
#**************************************************************************************
.PHONY: samtools picard-tools bedtools

seqtk:
	git clone https://github.com/lh3/seqtk.git
	cd seqtk && $(MAKE_NOFLAGS)
	$(call LINK_TO_BIN,`pwd`/seqtk/seqtk)

tabtk:
	git clone https://github.com/lh3/tabtk.git
	cd tabtk && $(MAKE_NOFLAGS)
	$(call LINK_TO_BIN,`pwd`/tabtk/tabtk)

samtools: samtools/1.2
samtools/1.2: VERSION=1.2
samtools/1.2:
	wget -N https://github.com/samtools/samtools/releases/download/$(VERSION)/samtools-$(VERSION).tar.bz2
	tar -xjf samtools-*.tar.bz2
	mkdir -p samtools && mv samtools-*/ $@
	cd $@ && make

picard-tools: picard-tools/1.135
picard-tools/1.135: VERSION=1.135
picard-tools/1.135:
	wget -N https://github.com/broadinstitute/picard/releases/download/$(VERSION)/picard-tools-$(VERSION).zip
	unzip picard-tools-*.zip
	mkdir picard-tools && mv picard-tools-*/ $@
	rm picard-tools-*.zip

bedtools: bedtools/2.24.0
bedtools/2.24.0: VERSION=2.24.0
bedtools/2.24.0:
	wget -N https://github.com/arq5x/bedtools2/releases/download/v$(VERSION)/bedtools-$(VERSION).tar.gz
	tar -xzf bedtools-*.tar.gz
	mkdir -p bedtools && mv bedtools2 $@
	cd $@ && $(MAKE_NOFLAGS)
	rm bedtools-*.tar.gz

#**************************************************************************************
#******************        QC and QF       ********************************************
#**************************************************************************************
.PHONY: FastQC jellyfish jellyfish2 FLASH fqtrim

FastQC: FastQC/0.11.3
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

jellyfish: jellyfish/1.1.11
jellyfish/1.1.11: VERSION=1.1.11
jellyfish/1.1.11: 
	-wget -N http://www.cbcb.umd.edu/software/jellyfish/jellyfish-$(VERSION).tar.gz	
	tar -xzf jellyfish-$(VERSION).tar.gz
	mkdir -p jellyfish 
	cd jellyfish-*/ && ./configure --prefix=$(shell pwd)/$@ && make && make install
	rm -r jellyfish-$(VERSION) jellyfish-$(VERSION).tar.gz

jellyfish2: jellyfish2/2.2.3
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

fqtrim: fqtrim/0.94
fqtrim/0.94: VERSION=0.94
fqtrim/0.94:
	-wget -N http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-$(VERSION).tar.gz
	tar -xzf fqtrim-*.tar.gz
	mkdir -p fqtrim && mv fqtrim-*/ $@
	cd $@ && make release
	rm fqtrim-*.tar.gz

pandaseq:
	echo "Requires zlib, bzip2 and libtool"
	#sudo apt-get install zlib1g-dev libbz2-dev libltdl-dev libtool
	git clone http://github.com/neufeld/pandaseq.git/
	cd pandaseq && bash autogen.sh && ./configure --prefix=`pwd`/dist && make && make install && make clean
	$(call LINK_TO_BIN,`pwd`/pandaseq/dist/bin/*)

FLASH: FLASH/1.2.11
FLASH/1.2.11: VERSION:=1.2.11
FLASH/1.2.11:
	wget -N http://downloads.sourceforge.net/project/flashpage/FLASH-$(VERSION).tar.gz
	tar -xzf FLASH-*.tar.gz
	mkdir -p FLASH && mv FLASH-*/ $@
	cd $@ && make
	rm FLASH-*.tar.gz

#**************************************************************************************
#******************        READ MAPPERS/ALIGNERS     **********************************
#**************************************************************************************
.PHONY: bwa bowtie2 last

bwa: bwa/0.7.12
bwa/0.7.12: VERSION=0.7.12
bwa/0.7.12:
	wget -N http://downloads.sourceforge.net/project/bio-bwa/bwa-$(VERSION).tar.bz2
	tar -xjf bwa-*.tar.bz2
	mkdir -p bwa
	mv bwa-*/ $@
	cd $@ && make
	rm bwa-*.tar.bz2

bowtie2: bowtie2/2.2.5
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

last: last/572
last/572: VERSION=572
last/572:
	wget -N http://last.cbrc.jp/last-$(VERSION).zip
	unzip last-$(VERSION).zip
	mkdir -p $@
	cd last-* && make && make install prefix=../$@
	rm -r last-*.zip  last-*/
	 
#**************************************************************************************
#******************        ASSEMBLERS     *********************************************
#**************************************************************************************
.PHONY: fermi SPAdes megahit kmc

fermi: fermi/1.1
fermi/1.1: VERSION=1.1
fermi/1.1:
	wget -N https://github.com/downloads/lh3/fermi/fermi-$(VERSION).tar.bz2
	tar -xjf fermi-*
	mkdir -p fermi && mv fermi-*/ $@
	cd $@ && make
	rm fermi-*.tar.bz2

megahit: megahit/0.3.3
megahit/0.3.3: VERSION=0.3.3
megahit/0.3.3:
	wget -N https://github.com/voutcn/megahit/releases/download/v$(VERSION)/megahit_v$(VERSION)_LINUX_CUDA6.5_sm350_x86_64-bin.tar.gz
	tar -xzf megahit_v*.tar.gz
	mkdir -p megahit && mv megahit_v*/ $@
	rm megahit_*.tar.gz

SPAdes: SPAdes/3.5.0
SPAdes/3.5.0: VERSION=3.5.0
SPAdes/3.5.0:
	wget -N http://spades.bioinf.spbau.ru/release$(VERSION)/SPAdes-$(VERSION)-Linux.tar.gz
	tar -xzf SPAdes*.tar.gz
	mkdir -p SPAdes && mv SPAdes-*/ $@
	rm SPAdes-*-Linux.tar.gz

#Assumes python3 installed from anaconda in a virtualenv called py3k
Fastaq:
	git clone https://github.com/sanger-pathogens/Fastaq.git
	cd Fastaq && $(python3_bin) setup.py install

kmc: kmc/2.2.0
kmc/2.2.0: VERSION:=2.2.0
kmc/2.2.0:
	mkdir -p $@
	cd $@ && wget -N http://sun.aei.polsl.pl/REFRESH/kmc/downloads/2.2.0/linux/kmc
	cd $@ && wget -N http://sun.aei.polsl.pl/REFRESH/kmc/downloads/2.2.0/linux/kmc_dump

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

hmmer3/3.1b1: VERSION=3.1b1
hmmer3/3.1b1:
	wget -N ftp://selab.janelia.org/pub/software/hmmer3/$(VERSION)/hmmer-*-linux-intel-x86_64.tar.gz
	tar -xzf hmmer-*-linux-intel-*.tar.gz
	mkdir -p hmmer3/ && mv hmmer-*/ $@
	@echo "Binaries in $@/binaries"

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

metaphlan2/901cc5778eed: VERSION=901cc5778eed
metaphlan2/901cc5778eed:
	wget -N https://bitbucket.org/biobakery/metaphlan2/get/$(VERSION).zip
	unzip $(VERSION).zip
	mkdir -p metaphlan2 && mv biobakery-metaphlan2-$(VERSION) $@

#************************************************************************
#****************      VARIANT CALLING       ****************************
#************************************************************************
.PHONY: freebayes vcflib

freebayes: freebayes/0.9.21
freebayes/0.9.21: VERSION=0.9.21
freebayes/0.9.21:
	git clone --recursive git://github.com/ekg/freebayes.git
	cd freebayes && git checkout tags/v$(VERSION)
	mv freebayes $(VERSION) && mkdir -p freebayes && mv $(VERSION) freebayes
	cd $@ && $(MAKE_NOFLAGS)

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

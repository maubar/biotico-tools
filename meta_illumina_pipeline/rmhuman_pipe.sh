#!/bin/bash
#SBATCH -A b2011088

#SBATCH -p node -n 1

#SBATCH -t 4:00:00

#SBATCH -J rmhuman_d1
#SBATCH --mail-user mauricio.barrientos@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -e slurm-%j.err

module load bwa/0.7.4
module load samtools/0.1.19

READ_DIR=/proj/b2011088/viruswork/samples/illumina
OUT_DIR=/home/mauricio/glob/results/vm003/rmhuman

R1=$READ_DIR"/R1.fastq.gz"
R2=$READ_DIR"/R2.fastq.gz"
REF=/proj/b2011088/db/GRCh37/iGenome/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
OUT_PRE=$OUT_DIR"/P431_101"

#Output files
NESONI_PRE=$OUT_PRE"_q20h"
BWA_PRE=$NESONI_PRE"_grch37"
RMSEC_PRE=$BWA_PRE"_rmsec"
SORTSAM_PRE=$RMSEC_PRE"_sort"
FILTSAM_PRE=$SORTSAM_PRE"_unmapped"
SAM2FQ_PRE=$FILTSAM_PRE

DEBUG=false

#File checks!
if ! [ -e $R1 ]; then echo "No R1 file: $R1" >&1; exit 1; fi
if ! [ -e $R2 ]; then echo "No R2 file: $R2" >&1; exit 2; fi
if ! [ -d $OUT_DIR ]; then echo "No out dir: $OUT_DIR" >&1; exit 3; fi

#Quality filtering
time nesoni clip --homopolymers yes --quality 20 --length 75 --out-separate yes $NESONI_PRE pairs: $R1 $R2
echo $'\n***'  >&2

#Map to human
bwa mem -t 8 -T 20 -M $REF $NESONI_PRE"_R1.fq.gz" $NESONI_PRE"_R2.fq.gz" > $BWA_PRE"_pairs.sam"
echo $'\n***'  >&2
bwa mem -t 8 -T 20 -M $REF $NESONI_PRE"_single.fq.gz" > $BWA_PRE"_single.sam"
echo $'\n***'  >&2
if ! "$DEBUG" ; then rm "$NESONI_PRE"_R1.fq.gz "$NESONI_PRE"_R2.fq.gz "$NESONI_PRE"_single.fq.gz ; fi

#Remove bwa secondary alignments(if it maps to human with both pairs it is good to remove)
samtools view -F 256 -hSb -o $RMSEC_PRE"_pairs.bam" $BWA_PRE"_pairs.sam" 
echo $'\n***'  >&2
samtools view -F 256 -hSb -o $RMSEC_PRE"_single.bam" $BWA_PRE"_single.sam" 
echo $'\n***'  >&2
if ! "$DEBUG" ; then rm $BWA_PRE"_pairs.sam" $BWA_PRE"_single.sam"; fi

#Sort by queryname
run_picard SortSam.jar INPUT=$RMSEC_PRE"_pairs.bam" OUTPUT=$SORTSAM_PRE"_pairs.bam" SORT_ORDER=queryname
echo $'\n***'  >&2
run_picard SortSam.jar INPUT=$RMSEC_PRE"_single.bam" OUTPUT=$SORTSAM_PRE"_single.bam" SORT_ORDER=queryname
echo $'\n***'  >&2
if ! "$DEBUG" ; then rm $RMSEC_PRE"_pairs.bam" $RMSEC_PRE"_single.bam"; fi

#Get alignment flag stats
samtools view $SORTSAM_PRE"_pairs.bam" | cut -f 2,5  > $SORTSAM_PRE"_flag_mapq_pairs.txt"
echo $'\n***' >&2
samtools view $SORTSAM_PRE"_single.bam" | cut -f 2,5  > $SORTSAM_PRE"_flag_mapq_single.txt"
echo $'\n***' >&2


#Keep only reads that did not map confidently (with both pairs)
run_picard FilterSamReads.jar INPUT=$SORTSAM_PRE"_pairs.bam" OUTPUT=$FILTSAM_PRE"_pairs.bam" FILTER=excludeAligned SORT_ORDER=queryname WRITE_READS_FILES=False
echo $'\n***'  >&2
run_picard FilterSamReads.jar INPUT=$SORTSAM_PRE"_single.bam" OUTPUT=$FILTSAM_PRE"_single.bam" FILTER=excludeAligned SORT_ORDER=queryname WRITE_READS_FILES=False
echo $'\n***' >&2
if ! "$DEBUG" ; then rm $SORTSAM_PRE"_pairs.bam" $SORTSAM_PRE"_single.bam"; fi

#Convert unmapped reads to Fastq for assembly
run_picard SamToFastq.jar INPUT=$FILTSAM_PRE"_pairs.bam" FASTQ=$SAM2FQ_PRE"_R1.fq" SECOND_END_FASTQ=$SAM2FQ_PRE"_R2.fq"
echo $'\n***' >&2
run_picard SamToFastq.jar INPUT=$FILTSAM_PRE"_single.bam" FASTQ=$SAM2FQ_PRE"_single.fq" 
echo $'\n***' >&2
if ! "$DEBUG" ; then rm $FILTSAM_PRE"_pairs.bam" $FILTSAM_PRE"_single.bam"; fi

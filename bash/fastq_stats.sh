#!/bin/bash

#Outputs the read count and and base count of a fastq file
set -euo pipefail

FASTQ_FILE=$1
CAT_TOOL="cat"

if [ ${FASTQ_FILE: -3} == ".gz" ]; then CAT_TOOL="zcat"; fi

${CAT_TOOL} ${FASTQ_FILE} |  sed -n "2~4p" | awk 'BEGIN{bases=0;reads=0;} {bases+=length($1);reads+=1} END{print reads,bases}'

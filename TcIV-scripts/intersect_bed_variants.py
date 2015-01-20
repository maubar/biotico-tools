#/usr/bin/env python
#
#Intersects bedfiles with variant calls to look for variants common in both files
#
# Author: Mauricio Barrientos-Somarribas
# Email:  mauricio.barrientos@ki.se
#
# Copyright 2014 Mauricio Barrientos-Somarribas

import sys
import re

CHROM=0
START=1
STOP=2
SEQ_REF=5
SEQ_ALT=6
INFO=7
GT=9

def main():
	if len(sys.argv) == 3:
		bed_ref = open(sys.argv[1])
		bed_query = open(sys.argv[2])
	else:
		print("Program use:\n\t\tintersect_bed_variants.py <reference.bed> <query.bed>")
		raise Exception

	file_out = sys.stdout

	#Load reference variants
	ref_variants = {}
	for line in bed_ref:
		if line.lstrip("\n")[0] != "#":
			pos,genotype = extract_fields(line)
			ref_variants[pos] = genotype

	for line in bed_query:
		if line.lstrip("\n")[0] != "#":
			pos,genotype = extract_fields(line)
			if pos in ref_variants: #If position(chromosome, start,end) are equal
				ref_genotype = ref_variants[pos]
				if ref_genotype == genotype:
					file_out.write(line)

def extract_fields(line):
	fields = line.rstrip("\n").split("\t")
	pos = (fields[CHROM],fields[START],fields[STOP])
	gt = extract_GT(fields[GT])
	genotype = (fields[SEQ_REF],fields[SEQ_ALT],gt)
	return (pos,genotype)

def extract_GT(field_gt):
	comma_pos = field_gt.find(":")
	return field_gt[:comma_pos]

if __name__ == "__main__":
	main()

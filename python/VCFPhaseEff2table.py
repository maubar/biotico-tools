#1/usr/bin/env python

import sys
import re

CHROM=0
POS=1
INFO=7
GT=9

def main():
	if len(sys.argv) == 1:
		vcf_file = sys.stdin
	else:
		vcf_file = open(sys.argv[1])

	file_out = sys.stdout
	file_out.write("Chrom\tPos\tAF\tMQ\tGT\tEffect\tImpact\tGene_name\n")

	for line in vcf_file:
		if line.lstrip()[0] != "#":
			file_out.write(extract_fields(line))

def extract_fields(line):
	fields = line.rstrip("\n").split("\t")
	new_line = "\t".join( (fields[CHROM],fields[POS]))
	new_line += "\t"+ extract_AF(fields[INFO])
	new_line += "\t"+ extract_MQ(fields[INFO])
	new_line += "\t"+ extract_GT(fields[GT])
	new_line += "\t"+ "\t".join(extract_effect(fields[INFO]))
	


	return new_line+"\n"

def extract_GT(field_gt):
	comma_pos = field_gt.find(":")
	return field_gt[:comma_pos]


def extract_AF(field_info):
	af_pos = field_info.find("AF=")
	af_semicolon = field_info.find(";",af_pos)
	return field_info[af_pos+3:af_semicolon]

def extract_MQ(field_info):
	mq_pos = field_info.find("MQ=")
	mq_semicolon = field_info.find(";",mq_pos)
	return field_info[mq_pos+3:mq_semicolon]

def extract_effect(field_info):
	eff_pos = field_info.find("EFF=")
	result = ["-"] *3
	if eff_pos > -1:
		eff_semicolon = field_info.find(";",eff_pos)

		effect_line = field_info[eff_pos+4:eff_semicolon]
		eff_name_delim = effect_line.find("(")
		eff_name = effect_line[:eff_name_delim]
		eff_fields = effect_line[eff_name_delim+1:-1].split("|")	
		result = (eff_name,eff_fields[0],eff_fields[5])
		
	return result

if __name__ == "__main__":
	main()

#/usr/bin/env python
#
#Intersects the variants from VCF files based on position (CHROM, POS) and reports whether the genotypes match
#
# Author: Mauricio Barrientos-Somarribas
# Email:  mauricio.barrientos@ki.se
#
# Copyright 2015 Mauricio Barrientos-Somarribas


import sys
import re

CHROM = 0
POS = 1
REF = 3
ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT = 8
SAMPLE = 9


def main():
    if len(sys.argv) == 3:
        vcf_ref = open(sys.argv[1])
        vcf_query = open(sys.argv[2])
    else:
        print("Program use:\n\t\tintersect_vcf_variants.py <reference.vcf> <query.vcf>")
        raise Exception

    file_out = sys.stdout

    # Index variants from the query
    query_variants = {}
    for line in vcf_query:
        if line.rstrip("\n")[0] != "#":
            pos,genotype = extract_fields(line)
            query_variants[pos] = genotype

    # Report whether each variant of the reference is present or not in the sample
    file_out.write("CHR\tPOS\tQUERY_GT\tGT_MATCHES_REF\n")
    for line in vcf_ref:
        if line.rstrip("\n")[0] != "#":
            pos,genotype = extract_fields(line)
            if pos in query_variants:
                file_out.write("{}\t{}\t{}\t{}\n".format(pos[0],pos[1], query_variants[pos] ,  query_variants[pos] == genotype  ))


def extract_fields(line):
    fields = line.rstrip("\n").split("\t")
    pos = (fields[CHROM],int(fields[POS]))
    gt = extract_GT(fields[SAMPLE],fields[FORMAT])
    genotype = gt["GT"]
    return pos,genotype


def extract_GT(sample_col,format_col):
    gt = dict([(k,v) for k,v in zip(format_col.split(":"),sample_col.split(":"))])
    return gt

if __name__ == "__main__":
    try:
        main()
    except IOError as e:
        print e

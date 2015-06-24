#!/usr/bin/env python
"""
Script that takes a gff file and uses the output gff of the Interpro
annotation pipeline to substitute the transcript id for obtained
sequence names

Author: Mauricio Barrientos-Somarribas
Email:  mauricio.barrientos@ki.se

Copyright 2014 Mauricio Barrientos-Somarribas

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import sys
import argparse
import os.path

# Time of script execution
import time
import logging

import re

# *************CONSTANTS***********
# Interpro GFF columns
GENE_ID = 0
ANNOT_TOOL = 1
TAGS = 8

# **************GLOBAL VARIABLES ***********
putative_prot_id = 0

# ****************Begin of Main ***************
def main(argv):
    # 1) Load Interprot gene annotation
    interprot_annot = {}
    for line in argv.interprot_gff:
        if line[0] != "#": # ignore comment lines
            line_fields = line.rstrip("\n").split("\t")

            # Extract SMART or Pfam annotation, but prefer SMART to Pfam
            if line_fields[ANNOT_TOOL] in ("Pfam", "SMART"):
                # Extract gene id (Removing the .t1 suffix)
                dot_position = line_fields[GENE_ID].rfind(".")
                if dot_position != -1:
                    line_fields[GENE_ID] = line_fields[GENE_ID][:dot_position]
                # Prefer SMART over Pfam annotation
                if line_fields[ANNOT_TOOL] == "SMART" or (line_fields[GENE_ID] not in interprot_annot):
                    signature_desc = extract_signature_desc(line_fields[TAGS])
                    interprot_annot[line_fields[GENE_ID]] = signature_desc

        elif line.rstrip("\n") == "##FASTA":
            break

    # 2) Parse Augustus gtf file and annotate
    for line in argv.augustus_gtf:
        if line[0] == "#":
            argv.output_file.write(line)
        else:
            argv.output_file.write(process_gtf_line(line, interprot_annot))
# END OF MAIN

def extract_signature_desc(tag_line):
    signature_desc = ""
    tag_fields = tag_line.split(";")
    for field in tag_fields:
        if "signature_desc" in field:
            value_start = field.find("=")
            signature_desc = field[value_start+1:]
    return signature_desc

def get_signature_description(gene_id, signature_descriptions):
    global putative_prot_id
    return signature_descriptions[gene_id] if gene_id in signature_descriptions \
        else "putative_protein_{}".format(putative_prot_id)

# Substitutes the original ID, or parent ID for the signature_desc from interprot pfam if available
def process_gtf_line(line, signature_descriptions):
    global putative_prot_id
    new_line = line
    id_match = re.search(r'gene_id "(.+?)"(;|\n)', line)
    if id_match:
        gene_id = id_match.group(1)
        new_gene_id = get_signature_description(gene_id, signature_descriptions)
        if new_gene_id.startswith("putative_protein"):
            putative_prot_id += 1

        new_line = re.sub(r'gene_id "(.+?)"(;|\n)', 'gene_id="{}"\g<2>'.format(new_gene_id), new_line)
        new_line = re.sub(r'transcript_id "(.+?)"(;|\n)', 'transcript_id="{}"\g<2>'.format(new_gene_id), new_line)
    return new_line


def validate_args(_):
    return True


if __name__ == '__main__':
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Substitutes the gene_id and transcript_id attributes \
        of an augustus GTF output file with protein names using the gff output from the EBI Interpro pipeline")
    parser.add_argument("interprot_gff",help="Input file", type=file)
    parser.add_argument("augustus_gtf",help="Input file", type=file)

    parser.add_argument("-o", "--output-file", type=argparse.FileType("w"), default=sys.stdout, help="Name of the output file" )
    parser.add_argument("-l", "--log-file", default=None, help="Name of the log file")

    args = parser.parse_args()

    if validate_args(args):
        # Initialize log
        log_level = logging.INFO
        if args.log_file:
            logging.basicConfig(filename=args.log_file, level=log_level)
        else:
            logging.basicConfig(stream=sys.stderr, level=log_level)

        time_start = time.time()
        main(args)
        logging.info("Time elapsed: "+str(time.time() - time_start)+"\n")
    else:
        logging.error("Invalid arguments. Exiting script\n")
        sys.exit(1)
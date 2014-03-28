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

#Time of script execution
import time
import logging

from collections import OrderedDict,deque
import re

#*************CONSTANTS***********
#Interpro GFF columns
GENE_ID = 0
ANNOT_TOOL = 1
TAGS = 8

#**************GLOBAL VARIABLES ***********
putative_prot_id = 0

#****************Begin of Main ***************
def main(args):
	interprot_annot = {}
	for line in args.gff_interprot:
		if line[0] != "#" : #ignore comment lines
			line_fields = line.rstrip("\n").split("\t")
			
			#Extract from Pfam lines the
			if line_fields[ANNOT_TOOL].lower() == "pfam":
				signature_desc = extract_signature_desc(line_fields[TAGS])
				#Remove the .t1 to the gene_id
				dot_position = line_fields[GENE_ID].rfind(".")
				if dot_position != -1:
					line_fields[GENE_ID] = line_fields[GENE_ID][:dot_position]
				interprot_annot[  line_fields[GENE_ID] ] = signature_desc
		elif line.rstrip("\n") == "##FASTA":
			break

	for line in args.gff_to_annotate:
		if line[0] == "#":
			args.output_file.write(line)
		else:
			args.output_file.write(processAugLine(line, interprot_annot))
#END OF MAIN

def extract_signature_desc(tag_line):
	signature_desc = ""
	tag_fields = tag_line.split(";")
	for field in tag_fields:
		if "signature_desc" in field:
			value_start = field.find("=")
			signature_desc = field[value_start+1:]
	return signature_desc

def getSignatureDescription(gene_id, signature_descriptions):
	global putative_prot_id
	return signature_descriptions[gene_id] if gene_id in signature_descriptions else "putative_protein"+str(putative_prot_id)

#Substitutes the orignal ID, or parent ID for the signature_desc from interprot pfam if available
def processAugLine(line, signature_descriptions):
	global putative_prot_id
	new_line = line
	if "\tgene" not in line:
		id_match = re.search( r"ID=(.+?)(\..*?)?(;|$)", line)
		if id_match:
			gene_id = id_match.group(1)
			new_gene_id = getSignatureDescription(gene_id,signature_descriptions)
			new_line = re.sub(r"ID=(.+?)(\..+?)(;|\n)", "ID='"+new_gene_id+r"'\g<3>",line)
		if "\ttranscript" not in line:
			parent_match = re.search( r"Parent=(.+?)(\..*?)?(;|$)", line)
			if parent_match:
				parent_id = parent_match.group(1)
				new_parent_id = getSignatureDescription(parent_id,signature_descriptions)
				new_line = re.sub(r"Parent=(.+?)(;|\n)", "Parent='"+new_parent_id+r"'\g<2>",new_line)
	else:
		#New protein, increase counter
		putative_prot_id += 1
	return new_line


def validate_args(args):
	return True


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Substitutes the ID attributes of a gff file with \
		protein names using the gff output from the EBI Interpro pipeline")
	parser.add_argument("gff_interprot",help="Input file", type=file)
	parser.add_argument("gff_to_annotate",help="Input file", type=file)

	parser.add_argument("-o","--output-file", type=argparse.FileType("w"), default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")
			
	args = parser.parse_args()
	
	if validate_args(args):
		#Initialize log
		log_level = logging.INFO
		if args.log_file:
			logging.basicConfig(filename=args.log_file,level=log_level)
		else:
			logging.basicConfig(stream=sys.stderr,level=log_level)

		time_start = time.time()
		main( args )
		logging.info("Time elapsed: "+str(time.time() - time_start)+"\n")
	else:
		logging.error("Invalid arguments. Exiting script\n")
		sys.exit(1)


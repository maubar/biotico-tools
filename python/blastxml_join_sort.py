#!/usr/bin/env python
"""
Blast XML parser that extracts the best hits for each query sequence
and presents the results in order by e-value in
tab-separated format

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

#Time of script execution and logging module
import time
import logging
import re

import itertools
from collections import *

from Bio.Blast import NCBIXML

#****************Begin of Main ***************
def main(args):
	filenames = extract_blast_filenames( args.root_folder )
	#Filter only refseqvir
	filtered_filenames = filter(lambda x:"refseqvir" in x, filenames)
	#Where hits will be accumulated
	results = []
	for blast_file in filtered_filenames:
		logging.info("Processing file: "+blast_file)
		assembler, blast_type , database = extract_metadata_from_blast_filename(blast_file)
		with open(blast_file,"r") as blast_fh:
			iterations = NCBIXML.parse(blast_fh)

			logging.info("Parsing XML file...")
			for iteration_num,it in enumerate(iterations):
				for hit_num,hit in enumerate(it.alignments):
					hsp = hit.hsps[0] #Only first hsp
					#Store relevant fields in a tuple
					if hsp.expect < args.threshold:
						visual_alignment = None
						if args.alignments:
							visual_alignment = (hsp.sbjct_start, hsp.sbjct, hsp.sbjct_end, hsp.match, hsp.query_start,hsp.query,hsp.query_end )
						species = extract_species_name(hit.hit_def,blast_type)
						q_cov , pct_id, pct_pos = calc_hsp_stats(hsp, it.query_length)
						result = (hsp.expect,
								species,
								it.query_length,
								q_cov,
								pct_id,
								pct_pos,
								it.query,
								assembler,
								blast_type,
								database,
								visual_alignment)
						results.append(result)

					#Limit number of hits to report
					if hit_num >= (args.max_hits-1):
						break
				#Keep track of number of iterations processed
				if iteration_num % 2000 == 0 :
					logging.info("\tIterations processed: "+str(iteration_num)+"\r")
		logging.info("XML parsed! " + str(iteration_num) + " iterations\n")

	#Sort
	logging.info("Sorting hits by e-value...\n")
	results.sort(key=lambda x:x[0],reverse=False)
	logging.info("Sort finished!\n")

	#Write out sorted results to file
	output_columns = ["e-value", "species", "seq_len", "q_cov","pct_id","pct_pos","ctg_id", "asm", "tool", "database","other"]
	args.output_file.write("\t".join(output_columns)+"\n")
	for res in results:
		args.output_file.write( "\t".join([str(x) for x in res[:9]])+"\n")
		if args.alignments:
			#print alignments as well
				alignment = res[9]
				args.output_file.write( "\n"+formatAlignment(alignment)+"\n\n" )


#*****************End of Main**********************
#alignment receives a tuple with the following fields:
#	hsp.sbjct_start, hsp.sbjct, hsp.sbjct_end, hsp.match, hsp.query_start,hsp.query,hsp.query_end
def formatAlignment(alignment):
	ref_line = "\tRef: "+str(alignment[0])
	qry_line = "\tQry: "+str(alignment[4])
	separator = "\t"
	#Space padding
	ref_line += (" " * max(0,len(qry_line)-len(ref_line)-1) )+ separator
	qry_line += (" " * max(0,len(ref_line)-len(qry_line)-1) )+ separator

	midline = "\t"+(" " * (len(ref_line)-2) )+ separator

	#Add sequences , padding and end coordinates
	ref_line += alignment[1] + separator + str(alignment[2])
	midline += alignment[3]
	qry_line += alignment[5] + separator + str(alignment[6])

	return ref_line+"\n"+midline+"\n"+qry_line

def extract_species_name(hit_name,blast_type):
	species = hit_name
	if blast_type == 'blastn':
		species = hit_name.split(",")[0]
	elif blast_type in ['blastp','blastx']:
		species = hit_name.split("[")[1].rstrip("] \n")
	return species

#Calculates percent identity, percent positives and query coverage
def calc_hsp_stats( hsp , qseq_len):
	pct_id = 100.0 * hsp.identities / hsp.align_length
	pct_pos = 100.0 * hsp.positives / hsp.align_length
	qcov = 100.0 * (hsp.query_end - hsp.query_start) / qseq_len
	return pct_id, pct_pos, qcov


def extract_metadata_from_blast_filename(filename):
	basename = os.path.basename(filename)
	asm_idx = basename.find("_asm_")
	if asm_idx > -1:
		fields = basename[asm_idx+1:].split("_")
		assembler = fields[1]
		tool = fields[2]
		db = fields[3].split(".xml")[0]
#paired-end or singletons
	else:
		match = re.search(r"_([^_]+?)_(blast[nxp])_([^_]+?)\.xml$",basename)
		if match:
			assembler = match.group(1)
			tool = match.group(2)
			db = match.group(3)
	return assembler,tool,db

def extract_blast_filenames(root_path):
	filenames = []
	for root_folder,dirs,files in os.walk(root_path):
		if "blast" in root_folder:
			filenames += [os.path.join(root_folder,f) for f in files if f.endswith(".xml") ]
	return filenames

def validate_args(args):
	return True

if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="The program extracts the highest scoring hits for each query \
		sequence and outputs them ordered by e-value")

	parser.add_argument("root_folder",help="Folder in which to look for blast xml files", default="./")
	parser.add_argument("-o","--output-file", type=argparse.FileType('w'), default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")

	parser.add_argument("--max_hits",type=int,default="10",help="Max number of hits to report for each query sequence")
	parser.add_argument("-a","--alignments",action='store_true',help="Include a visualization for each alignment on the file")
	parser.add_argument("-t","--threshold",type=float,default=1,help="Exclude hits with e-value lower than threshold. Default is '1' (all hits)")

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

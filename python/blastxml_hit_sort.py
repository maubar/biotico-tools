#!/usr/bin/env python
"""
Blast XML parser that extracts the best hits for each query sequence
and presents the results in order by e-value or bitscore in
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

import itertools
from collections import *

from Bio.Blast import NCBIXML

#****************Begin of Main ***************
def main(args):
	results = []
	iterations = NCBIXML.parse(args.blast_xml_file)

	logging.info("Parsing XML file...\n")
	for it in iterations:
		hit_num = 0
		for hit in it.alignments:
			hsp = hit.hsps[0] #Only first hsp
			#Store relevant fields in a tuple
			if args.evalue:
				sorting_key = hsp.expect
				threshold = 0.001
			elif args.bitscore:
				sorting_key = hsp.bits

			if sorting_key < threshold:
				result = (sorting_key, it.query, it.query_length, hsp.align_length,hsp.identities,hsp.gaps, hit.hit_id, hit.hit_def, hsp.hseq , hsp.midline, hsp.qseq )
				results.append(result)
			#Limit number of hits to report
			hit_num +=1
			if hit_num == args.max_hits:
				break

	logging.info("XML parsed!\n")

	#Sort
	logging.info("Starting sort...\n")
	results.sort(key=lambda x:x[0],reverse=args.bitscore)

	logging.info("Sort finished!\n")

	sort_col =  "E-value" if args.evalue else "Bitscore"
	#Write out sorted results to file
	output_columns = [sort_col, "Query_id", "Query_length", "Alignment_length","Identities"," Gaps","Hit_ID", "Hit_name"]
	args.output_file.write("\t".join(output_columns)+"\n")
	for res in results:
		args.output_file.write( "\t".join([str(x) for x in res[:8]])+"\n\n" )
		args.output_file.write(res[8]+"\n")
		args.output_file.write(res[9]+"\n")
		args.output_file.write(res[10]+"\n\n")

#*****************End of Main**********************
def validate_args(args):
	return True

if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="The program extracts the highest scoring hits for each query \
		sequence and outputs them ordered by e-value")

	parser.add_argument("blast_xml_file",help="Blast XML file to parse",nargs="?", type=file, default=sys.stdin)
	parser.add_argument("-o","--output-file", type=file, default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")

	parser.add_argument("--max_hits",type=int,default="10",help="Max number of hits to report for each query sequence")

	sort_criteria = parser.add_mutually_exclusive_group(required=True)
	sort_criteria.add_argument("-b","--bitscore",action="store_true",help="Sort hits from highest to lowest bitscore")
	sort_criteria.add_argument("-e","--evalue",action="store_true", help="Sort hits from lowest to highest e-value")

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

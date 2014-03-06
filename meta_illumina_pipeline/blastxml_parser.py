#!/usr/bin/env python
"""
Program processes blast xml output, and outputs a summarized version
of the hit results , keeping only the relevant information for the
pipeline.


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

import re
import math

import itertools
from collections import * 

import pdb
import inspect


#****************Begin of Main ***************
def main(args):
	queries_processed = 0
	for line in args.blast_xml:
		if re.match(r"\s*<Iteration>\s*$",line):
			iteration_status = processIteration(args.blast_xml,args.output_file)
			assert(iteration_status)
			queries_processed +=1

	args.log_file.write("Total iterations processed: "+str(queries_processed)+"\n")


#*****************End of Main**********************	
def processIteration(blastxml_fh,output_fh):
	processing_succesful = False
	iteration_data = {}
	for line in blastxml_fh:
		#End of iteration, exit succesfully
		if "</Iteration>" in line:
			processing_succesful = True
			break
		#Start processing hits for query
		elif "<Iteration_hits>" in line:
			iteration_data["Hits"] = []
			hit_status = processHits(blastxml_fh,output_fh,iteration_data)
			#Check if there was a problem processing hits
			assert(hit_status)

		else: #Extract query attributes of interest
			pattern_match = re.search(r"<(Iteration_query-(def|len))>(?P<value>.+?)</\1>",line)
			if pattern_match:
				iteration_data["query-"+pattern_match.group(2)] = pattern_match.group("value")

	#Process blast hits
	#if "Hits" in iteration_data and iteration_data["Hits"]:
	output_fh.write(str(iteration_data)+"\n\n")
	return processing_succesful

def processHits(blastxml_fh, output_fh,iteration_data):
	processing_succesful = False
	for line in blastxml_fh:
		if "</Iteration_hits>" in line:
			processing_succesful = True
			break
		elif "<Hit>" in line:
			new_hit = {}
			hit_status = processHit(blastxml_fh,output_fh,new_hit)
			assert(hit_status)
			iteration_data["Hits"].append(new_hit)

	return processing_succesful

def processHit(blastxml_fh,output_fh,hit_data):
	processing_succesful = False
	for line in blastxml_fh:
		if "</Hit>" in line:
			processing_succesful = True
			break				
		elif "<Hit_hsps>" in line:
			hit_data["Hsps"] = []
			hsps_status = processHsps(blastxml_fh,output_fh, hit_data )
			assert(hsps_status)
		else: #Extract hit attributes
			hit_attr_match = re.search(r"<Hit_(.+?)>(.+?)</Hit_\1>",line)
			if hit_attr_match:
				hit_data[ hit_attr_match.group(1) ] = hit_attr_match.group(2)
	return processing_succesful

def processHsps(blastxml_fh,output_fh, hit_data):
	processing_succesful = False
	for line in blastxml_fh:
		if "</Hit_hsps>" in line:
			processing_succesful = True
			break
		elif "<Hsp>" in line:
			new_hsp = {}
			hsp_status = processHsp(blastxml_fh,output_fh,new_hsp)
			assert(hsp_status)
			hit_data["Hsps"].append(new_hsp)
	return processing_succesful

def processHsp(blastxml_fh,output_fh,hsp_data):
	processing_succesful = False
	for line in blastxml_fh:
		if "</Hsp>" in line:
			processing_succesful = True
			break	
		else: #Extract hsp attributes
			hsp_attr_match = re.search(r"<Hsp_(.+?)>(.+?)</Hsp_\1>",line)
			if hsp_attr_match:
				hsp_data[ hsp_attr_match.group(1) ] = hsp_attr_match.group(2)
	return processing_succesful

def validate_args(args):
	return True


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="WRITE DESCRIPTION HERE",epilog= "TEXT AFTER THE HELP")
	parser.add_argument("blast_xml",help="Input file",nargs="?", type=file, default=sys.stdin)
	
	parser.add_argument("-o","--output-file", type=file, default=sys.stdout, help="Name of the output file" )
	#parser.add_argument("-o","--output-prefix", default=None, help="Prefix of the output file(s)" )
	#parser.add_argument("-o","--output-folder", default="./", help="Folder where to output the results" )

	parser.add_argument("-l","--log-file",type=argparse.FileType("w"), default=sys.stderr, help="Name of the log file")
	parser.add_argument("-v","--verbose", action="store_true", help="Print extra information about the execution")
			
	args = parser.parse_args()
	
	if validate_args(args):
		time_start = time.time()
		try:
			main( args )
		except IOError as ioerr:
			if ioerr.errno != 32: #Broken pipe error(when using head)
				args.log_file.write(str(ioerr)+"\n")
				raise
			#args.log_file.write(str(ioerr)+"\n")
		except KeyboardInterrupt:
			args.log_file.write("Execution halted by user")
		else:
			elapsed_time = int ((time.time() - time_start)*1000) * 0.001
			args.log_file.write("Time elapsed: "+str(elapsed_time)+" seconds\n")
	else:
		sys.stderr.write("Invalid arguments. Exiting script\n")
		sys.exit(1)
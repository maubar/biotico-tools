#!/usr/bin/env python
"""
Download fastas from a list of gi identifiers of the form
gi|xxxx| from a specified NCBI database


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

import requests
import argparse
import os.path
import sys

#Time of script execution and logging module
import time
import logging

import re
#****************Begin of Main ***************
def main(args):
	#Read gis from file
	gi_regex = re.compile(r"gi\|(.+?)\|")
	numeric_gi_regex = re.compile(r"^[0-9]+$")

	processed_gi_list = []
	for line in args.input_file:
		#If id is just a number, assume it is a gi
		if numeric_gi_regex.match(line):
			processed_gi_list.append(line.rstrip("\n"))
		else: #If id has a gi|xxx| structure, extract the gi number
			regex_hit = gi_regex.search(line)
			if regex_hit:
				processed_gi_list.append(regex_hit.group(1))

	#Send the request to the NCBI via eutils
	gis_to_get = ",".join(set(processed_gi_list))
	logging.info("Extracting [{}] from NCBI {}".format(gis_to_get,args.db))
	r = requests.get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&id={}&rettype=fasta&retmode=text".format(args.db,gis_to_get))
	if r.status_code == 200:
		try:
			args.output_file.write(str(r.content.decode()))
		except Exception as e:
			logging.error("There was an error writing the output")
			logging.error(e.str())
	else:
		logging.error("The request was returned with status code {}".format(r.statuscode))

#*****************End of Main*********************

def validate_args(args):
	return True


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="WRITE DESCRIPTION HERE",epilog= "TEXT AFTER THE HELP")
	parser.add_argument("input_file",help="Input file",nargs="?", type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument("-o","--output-file", type=argparse.FileType('w'), default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")
	parser.add_argument("--db", default="nuccore", help="Name of the log file")
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

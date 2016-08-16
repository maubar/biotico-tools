#!/usr/bin/env python
import sys
import argparse
import os.path
import time
import logging

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def main(args):
	with open(args.fasta_input,"r") as input_f:
		total_reads = 0
		match_count = 0

		seq_id_set  = createSeqIdSet(args.sequence_ids)

		include_sequence = False
		for line in input_f:
			if is_header(line):
				include_sequence = (extractReadName(line) in seq_id_set) ^ args.exclude_seqs
				total_reads += 1
				if include_sequence:
					match_count += 1
					args.output_file.write(line)
			else:
				if include_sequence:
					args.output_file.write(line)

		logging.info("")
		logging.info( "****Stats****")
		logging.info( "Size of filtering set:\t{}".format( len(seq_id_set)) )
		logging.info( "Total reads:\t{}".format(total_reads) )
		logging.info( "Reads kept:\t{}".format(match_count) )
		logging.info( "Reads discarded:\t{}".format( total_reads - match_count) )

#*****************End of main**********************
def is_header(line):
	return line[0] ==  ">"

def extractReadName(line):
	return line[1:]#.rstrip("\n")

#******************I/O functions***********************************

def createSeqIdSet(seq_id_file):
	return set( seq_id_file.readlines() );

#*****************Bureocracy functions**************

#Drops path and extension(if exists
def extract_filename(filepath):
	filename = filepath.rstrip(" ").split("/")[-1]
	extension_begin = filename.rfind(".")
	if extension_begin > 0:
		filename_noext = filename[:extension_begin]
	else:
		filename_noext = filename
	return filename_noext

def validate_args(args):
	if args.fasta_input and not os.path.isfile(args.fasta_input):
		sys.stderr.write("Error! "+args.fasta_input+" does not exist!\n")
		sys.exit(1)

	return True

if __name__ == '__main__':
	log_level = logging.INFO
	logging.basicConfig(stream=sys.stderr,level=log_level)

	#Process command line arguments
	parser = argparse.ArgumentParser(description="Filters a fasta file based on the specified sequence_ids file")
	parser.add_argument("fasta_input",help="Fasta file to filter")
	parser.add_argument("sequence_ids",type=argparse.FileType("r"), help="List of sequences to whitelist/blacklist")
	parser.add_argument("-o","--output-file", type=argparse.FileType("w"), default=sys.stdout, help="Name of the output file. Defaults to stdout" )
	parser.add_argument("-e","--exclude-seqs",default=False,action="store_true",help="Treats the sequence_ids as a blacklist instead of a whitelist")

	args = parser.parse_args()

	if validate_args(args):
		time_start = time.time()
		main( args )
		logging.info("Time elapsed: "+str(time.time() - time_start)+"\n")

	else:
		logging.error("Invalid arguments. Exiting script\n")
		sys.exit(1)

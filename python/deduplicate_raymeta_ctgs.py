#!/usr/bin/env python
import sys
import argparse
import os.path
import time

def main(args):
	with open(args.fasta_input,"r") as input_f:
		seq_id_set  = set()
		is_duplicated = False
		unique_seqs = 0
		for line in input_f:
			if is_header(line):
				seq_id = extractReadName(line)
				is_duplicated = seq_id in seq_id_set
				if not is_duplicated:
					unique_seqs += 1
					seq_id_set.add( seq_id )
					args.output_file.write(line)
				else:
					#If one sequence repeats, all the rest will be duplicates
					break
			else:
				if not is_duplicated:
					args.output_file.write(line)
		args.output_file.close()
		sys.stderr.write( "****Stats****\n ")
		sys.stderr.write( "Kept reads:\t"+str(unique_seqs)+"\n" )
		

#*****************End of main**********************
def is_header(line):
	return line[0] ==  ">"

def extractReadName(line):
	return line[1:]#.rstrip("\n")

#*****************Bureocracy functions**************

#Drops path and extension(if exists)
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

	if not args.output_file:
		args.output_file = sys.stdout
	else:
		args.output_file = open( args.output_file, "w")

	return True

if __name__ == '__main__':
	start_time = time.time()

	#Process command line arguments
	parser = argparse.ArgumentParser(description="Filters a fasta file based on the specified headers on a file")
	parser.add_argument("fasta_input",help="Fasta file to filter")
	parser.add_argument("-o","--output-file",help="Name for the output file. Default: stdout")

	args = parser.parse_args()

	if validate_args(args):
		main( args )
		sys.stderr.write(str(time.time() - start_time) + "seconds\n")

	else:
		sys.stderr.write("Invalid arguments. Exiting script")
		sys.exit(1)

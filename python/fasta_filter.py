#!/usr/bin/env python
import sys
import argparse
import os.path
import time
#from HTSeq import GenomicInterval 
#from collections import defaultdict

def main(args):
	with open(args.fasta_input,"r") as input_f:
		output_file = createOutputFile(args.output_file,args.exclude_seqs)
		
		total_reads = 0
		match_count = 0
		
		header_set  = createHeaderSet(args.header_file)
		
		include_sequence = False
		for line in input_f:
			if is_header(line): 
				include_sequence = (extractReadName(line) in header_set) ^ args.exclude_seqs
				total_reads += 1
				if include_sequence:
					match_count += 1
					output_file.write(line)
			else:				
				if include_sequence:
					output_file.write(line)
							
		print
		print "****Stats****"
		print "Size of filtering set: "+str(len(header_set))
		
		print "Total reads:\t"+str(total_reads)
		print "Reads kept:\t"+str(match_count)
		print "Reads filtered out:\t"+str(total_reads - match_count)
		
#*****************End of main**********************	
def is_header(line):
	return line[0] ==  ">"

def extractReadName(line):
	return line[1:]#.rstrip("\n")

#******************I/O functions***********************************
def createOutputFile(prefix,exclude):
	output_file = None
	try:
		if exclude:
			output_file = open(prefix+"_filt.fa","w")
		else:
			output_file = open(prefix+"_excl.fa","w")
	
	except Exception as e:
		sys.stderr.write("An error ocurred creating the output file\n\n" + str(e))
		output_file = None
	return output_file

def createHeaderSet(header_file):
	return set( header_file.readlines() );

#def writeLine(line,file_handle):
#	file_handle.write(line)


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
	
	#if args.header_file and not os.path.isfile(args.header_file):
	#	sys.stderr.write("Error! "+args.header_file+" does not exist!\n")
	#	sys.exit(1)

	if not args.output_file:
		args.output_file = extract_filename(args.fasta_input)
		
	return True

if __name__ == '__main__':
	start_time = time.time()
	
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Filters a fasta file based on the specified headers on a file")
	parser.add_argument("fasta_input",help="Fasta file to filter")
	parser.add_argument("header_file",type=argparse.FileType("r"), help="Blast output file to process in tabular format")
	parser.add_argument("-o","--output-file",help="Name for the output file. Default: 'fasta_input' + _filtered ")
	parser.add_argument("-e","--exclude-seqs",action="store_true",help="Changes the default behavior of the program. \
								Instead of keeping the headers from the file, keeps what does not match.")
	
	args = parser.parse_args()
	
	if validate_args(args):
		main( args )
		print time.time() - start_time, "seconds"
		
	else:
		sys.stderr.write("Invalid arguments. Exiting script")
		sys.exit(1)

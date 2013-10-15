#!/usr/bin/env python
import sys
import argparse
import os.path
import time

def main(args):
	with open(args.input_reads,"r") as input_f:
		
		output_files = createOutputFiles(args.output_folder+args.output_prefix)
		
		total_read_count = 0
		single_read_count = 0
		
		expecting_pair = False
		first_pair = None
		
		current_read = getReadData(input_f)
		while( current_read ):
			if(not expecting_pair ):
				if(current_read["type"] == "S"):
					writeRead(current_read, output_files["S"] )
					single_read_count += 1
					total_read_count +=1
				else:
					expecting_pair = True
					first_pair = current_read
				
			else: #Previous was rev or Single
				if(current_read["type"] == "S" or first_pair["type"] == current_read["type"] ):
					error_read_line = 4*(2*total_read_count - single_read_count ) +1
					sys.stderr.write("Error, an unpaired read was found on line"+str(error_read_line)+"\n\n")
					sys.stderr.write("Previous read: "+str(first_pair)+"\n")
					sys.stderr.write("Current read: "+str(current_read)+"\n")
					raise Exception 
				else:
					writeRead( first_pair  , output_files[ first_pair["type"] ])
					writeRead( current_read, output_files[ current_read["type"] ])
					total_read_count += 1
					expecting_pair = False
					first_pair = None

			current_read = getReadData(input_f)

		print "***Stats***"
		print "Total read count "+str(total_read_count)
		print "Single read counts "+str(single_read_count)
		print "Paired-end read counts "+str(total_read_count - single_read_count)
		
#*****************End of main**********************	

#******************I/O functions***********************************
def createOutputFiles(prefix):
	output_files = {}
	try:
		output_files["S"] = open(prefix+"_single.fastq","w")
		output_files["F"] = open(prefix+"_R1.fastq","w")
		output_files["R"] = open(prefix+"_R2.fastq","w")
	except Exception as e:
		sys.stderr.write("An error ocurred creating the output files\n" + str(e))
		output_files = None
	return output_files

def writeRead(fq_read,file_handle):
	file_handle.write(fq_read["raw"])

def getReadData(file_handle):
	fq_read = {}
	read_header = file_handle.readline()
	if(read_header):
		#Process header
		header_fields = read_header.rstrip("\n").split(" ")
		fq_read["id"]= header_fields[0]
		for kv_pair in header_fields[1:]:
			key,value = kv_pair.split("=")
			fq_read[key] = value
		#Store raw read
		fq_read["raw"] = read_header+file_handle.readline()+file_handle.readline()+file_handle.readline()
		
		if("dir" in fq_read):
			fq_read["type"] = fq_read["dir"]
		else:
			fq_read["type"] = "S"
	else:
		fq_read = None
	
	return fq_read

def extract_filename(filepath):
	filename = filepath.rstrip(" ").split("/")[-1]
	extension_begin = filename.rfind(".")
	if extension_begin > 0 and filename[extension_begin:] in [".fastq",".fq"]:
		filename_noext = filename[:extension_begin]
	else:
		sys.stderr.write( "Unrecognized input file extension")
		raise Exception
	return filename_noext
	
	

def validate_args(args):

	if args.input_reads and not os.path.isfile(args.input_reads):
		sys.stderr.write("Error! "+args.input_reads+" does not exist!\n")
		sys.exit(1)
		
	if args.output_folder:
		if not os.path.isdir(args.output_folder):
			sys.stderr.write("Error! "+args.output_folder+" does not exist!\n")
			sys.exit(2)		
		if args.output_folder[-1] != "/":
			args.output_folder = args.output_folder + "/"
	else:
		args.output_folder = "./"
		
	if not args.output_prefix:
		args.output_prefix = extract_filename(args.input_reads)
		
	return True

if __name__ == '__main__':
	start_time = time.time()
	
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Splits the reads of an aggregated 454 fastq file into single-ends and paired-ends")
	parser.add_argument("input_reads",help="Fastq file to split")
	parser.add_argument("-o","--output-folder",help="Define folder for output files. Default is ./")
	parser.add_argument("-p","--output-prefix",help="Change name prefix for output files. Default is input_reads file name")
	args = parser.parse_args()
	
	if validate_args(args):
		main( args )
		print time.time() - start_time, "seconds"
		
	else:
		sys.stderr.write("Invalid arguments. Exiting script")
		sys.exit(1)

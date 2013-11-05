#!/usr/bin/env python
import sys
import argparse
import os.path
import time

def main(args):
	with open(args.input_reads,"r") as input_f:

		output_files = createOutputFiles(args.output_folder+args.output_prefix)

		total_read_count = 0

		expecting_pair = False
		first_pair = None

		current_read = getReadData(input_f)
		while( current_read ):
			if not expecting_pair:
				expecting_pair = True
				first_pair = current_read

			else:
				if first_pair["pair"] != current_read["pair"] :
					writeRead( first_pair  , output_files[ first_pair["pair"] ])
					writeRead( current_read, output_files[ current_read["pair"] ])
					total_read_count += 1
					expecting_pair = False
					first_pair = None
				else: #If both are R1 or R2 , ERROR!
					error_read_line = 4*(2*total_read_count - single_read_count ) +1
					sys.stderr.write("Error, an unpaired read was found on line"+str(error_read_line)+"\n\n")
					sys.stderr.write("Previous read: "+ str(first_pair)  +"\n")
					sys.stderr.write("Current read: " + str(current_read)+"\n")
					raise Exception

			current_read = getReadData(input_f)

		print "***Stats***"
		print "Total pair count "+str(total_read_count)

#*****************End of main**********************

#******************I/O functions***********************************
def createOutputFiles(prefix):
	output_files = {}
	try:
		output_files["1"] = open(prefix+"_R1.fastq","w")
		output_files["2"] = open(prefix+"_R2.fastq","w")
	except Exception as e:
		sys.stderr.write("An error ocurred creating the output files\n" + str(e))
		output_files = None
	return output_files

def writeRead(fq_read,file_handle):
	file_handle.write("@"+fq_read["id"]+"_"+fq_read["ref"]+
							" "+(":".join(map(str,fq_read["pos"])))+
							"_"+(":".join(map(str,fq_read["errors"])))+
							"/"+fq_read["pair"]+"\n"
							)
	file_handle.write(fq_read["seq"]+"\n+\n"+fq_read["qual"]+"\n")

def getReadData(file_handle):
	fq_read = {}
	read_header = file_handle.readline().rstrip("\n")
	if(read_header):
		#ASSERT
		if read_header[0]!="@":
			raise Exception("Unknown header")

		header_fields = read_header[1:].split(" ") #[1:] to remove the @

		if header_fields[0].find("/") == -1:
			raise Exception("Unable to find forward/reverse pair token")

		fq_read["id"],fq_read["pair"] = header_fields[0].split("/")

		if "reference=" not in header_fields[1]:
			raise Exception("No reference present in header")

		fq_read["ref"] = header_fields[1][10:]

		if "position=" not in header_fields[2]:
			raise Exception("No position in header")

		fq_read["pos"] = extract_pos(header_fields[2][9:])

		if len(header_fields) >= 4 and header_fields[3] !="":
			if "errors=" not in header_fields[3] and "description" not in header_fields[3]:
				print header_fields
				print read_header
				raise Exception("No errors field in header")
			fq_read["errors"] = count_errors( header_fields[3][7:])
		else:
			fq_read["errors"] = (0,0,0)

		#Store raw read
		fq_read["seq"] = file_handle.readline().rstrip("\n")

		if file_handle.readline().rstrip("\n") != "+":
			raise Exception("Malformed fastq file")

		fq_read["qual"] = file_handle.readline().rstrip("\n")

	else:
		fq_read = None

	return fq_read

def extract_pos(pos_string):
	strand = "fw"
	if "complement" in pos_string:
		strand = "rv"
		pos_string = pos_string[11:-1] #Remove 'complement(' and ')' at the end

	begin, end = pos_string.split("..")

	return (begin,end, strand)

def count_errors(error_string):
	#Substitutions , insertions, deletions
	return( error_string.count("%") , error_string.count("+"), error_string.count("-") )

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
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=
  """
  Parses the interleaved paired-end read output of Grinder to produce individual files for
  each paired-end and a less verbose header.

  The current header format is composed of:

  @<id>_<seq-id>_<start-pos>:<end-pos>:<strand>_<substitution-count>:<insertion-count>:<deletion-count>/<paired-end>
  """)
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

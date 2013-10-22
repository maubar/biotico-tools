#!/usr/bin/env python
import sys
import argparse
import os.path
import time
from HTSeq import GenomicInterval 
#from collections import defaultdict

def main(args):
	with open(args.input_file,"r") as input_f:
		
		output_files = createOutputFiles(args.output_folder+args.output_prefix,args.match_only,args.no_match_only)
		
		match_count = 0
		query_count = 0
		
		query_line = getQueryLine(input_f)
		current_qname = query_line["qname"]
		current_q_hits = {}
		
		while( query_line ):

			if query_line["qname"] != current_qname:
				query_count += 1
			
				#Define if the hits are considered a match or a no-match
				if hits_are_significant(current_q_hits):
					match_count += 1
					if not args.no_match_only:
						#Write to match
						writeToFile(current_qname, output_files["match"])
						#writeHitLog(current_qname, current_q_hits, "
						
				elif not args.match_only:
					#Write to no-match
					writeToFile(current_qname, output_files["nomatch"])
		  	
				#Update currents
				current_qname = query_line["qname"]
				current_q_hits = {}
			
			#Add current query line temporary list
			addQueryHit( query_line, current_q_hits)
			query_line = getQueryLine(input_f)

		#Process last query out of the while loop
		query_count += 1
			
		#Define if the hits are considered a match or a no-match
		if hits_are_significant(current_q_hits):
			match_count += 1
			if not args.no_match_only:
				#Write to match
				writeToFile(current_qname, output_files["match"])

		elif not args.match_only:
			#Write to no-match
			writeToFile(current_qname, output_files["nomatch"])
		
		print "***Stats***"
		print "Total queries:\t"+str(query_count)
		print "Queries with a valid hit:\t"+str(match_count)
		print "Queries with no valid hits:\t"+str(query_count - match_count)
		
#*****************End of main**********************	

def hits_are_significant(current_hits):
	for hit in current_hits:
		if type(current_hits[hit]) != list: #it was a length 50+ match, 
			return True
		elif len(current_hits[hit]) >= 2 : #if there are at least to distinct hits for a same protein
			#print current_hits[hit]
			return True
	
	return False

def addQueryHit(query_line,current_hits):
	hit_coords = query_line["subject_coord"]
	
	if query_line["subject"] not in current_hits:
		current_hits[query_line["subject"]] = []
	
	known_hits = current_hits[query_line["subject"]]
	
	if type(known_hits) == list:
		#ctal defined thresholds(arbitrary)
		if query_line["subject_coord"].length >= 50 and query_line["pct_id"] >= 45.0:   
			current_hits[query_line["subject"]] = query_line["subject_coord"]
		
		elif query_line["subject_coord"].length < 50 and query_line["pct_id"] >= 80.0:
			is_contained = False
			hits_to_pop = []
			for idx,known_hit in enumerate(known_hits):
				if query_line["subject_coord"].contains( known_hit):
					print str(query_line["subject_coord"])+" contains "+str(known_hit)
					hits_to_pop.append(idx)
				if query_line["subject_coord"].is_contained_in( known_hit):
					print str(query_line["subject_coord"])+" is contained by "+str(known_hit)
					is_contained = True
					break
			#Possible problem!
			if is_contained and len(hits_to_pop) != 0: 
				print "Possible problem !"
				
			for contained_hit in sorted(hits_to_pop,reverse=True):
				known_hits.pop(contained_hit)
				
			if not is_contained:
				known_hits.append( query_line["subject_coord"] )
				
		#else:
			#print "Rejected based on length and quality: len "+str(query_line["subject_coord"].length)+";  % id "+str(query_line["pct_id"])


#******************I/O functions***********************************
def createOutputFiles(prefix, match_only=False, no_match_only=False):
	output_files = {}
	try:
		if not no_match_only:
			output_files["match"] = open(prefix+"_match.txt","w")
		
		if not match_only:
			output_files["nomatch"] = open(prefix+"_nomatch.txt","w")
	
	except Exception as e:
		sys.stderr.write("An error ocurred creating the output files\n" + str(e))
		output_files = None
	return output_files

def writeToFile(query_name,file_handle):
	file_handle.write(query_name+"\n")

def writeHitLog(query_hits,file_handle):
	pass

def getQueryLine(file_handle):
	query_data = None
	query_line = file_handle.readline()
	if(query_line):
		query_data = {}
		#Process header
		query_fields = query_line.rstrip("\n").split("\t")
		
		query_data["qname"] = query_fields[0]
		query_data["subject"] = query_fields[1]
		query_data["pct_id"] = float(query_fields[2])
		#query_data["ali_len"] = int(query_fields[3])
		query_data["subject_coord"] = GenomicInterval(query_fields[1],int(query_fields[8]), int(query_fields[9])+1 , "+" )
		#Store raw line
		#query_data["raw_line"] = query_fields
		
	
	return query_data

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
	if args.input_file and not os.path.isfile(args.input_file):
		sys.stderr.write("Error! "+args.input_file+" does not exist!\n")
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
		args.output_prefix = extract_filename(args.input_file)
		
	return True

if __name__ == '__main__':
	start_time = time.time()
	
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Splits the output of blast into two disjunctive sets: queries with a 'good enough' hit as defined by the parameters, and queries with a poor or no hit")
	parser.add_argument("input_file",help="Blast output file to process in tabular format")
	parser.add_argument("-o","--output-folder",help="Define a folder for the output files. Default: './'")
	parser.add_argument("-p","--output-prefix",help="Change name prefix for output files. Default is input_reads file name")
	
	g = parser.add_mutually_exclusive_group()
	g.add_argument("-M", "--match-only",action="store_true",help="Only output the list of proteins with good hits")
	g.add_argument("-N", "--no-match-only", action="store_true",help="Only output the list of proteins without a proper match")
	
	args = parser.parse_args()
	
	if validate_args(args):
		main( args )
		print time.time() - start_time, "seconds"
		
	else:
		sys.stderr.write("Invalid arguments. Exiting script")
		sys.exit(1)

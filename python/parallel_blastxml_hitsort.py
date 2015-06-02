#!/usr/bin/env python
import sys
import argparse
import os.path

#Time of script execution and logging module
# import time
# import logging
import re
import itertools
# from collections import *
# import random

import multiprocessing as mp
from Bio.Blast import NCBIXML

def main():
	args = create_argparser()
	#Get files to process
	files_to_process = get_files_to_process(args)
	sys.stderr.write("Files to process:\n\t{}\n".format("\n\t".join(files_to_process)))

	#Get subset of sequences to extract hits from
	sequence_subset = load_sequence_subset(args.subset_file) if args.subset_file else None

	#Run blastxml parser and extract data
	process_pool = mp.Pool(args.num_proc)
	hits_per_file = process_pool.map(blastxml_parser(args,sequence_subset), files_to_process)
	#Merge/sort hits by e-value
	sys.stderr.write("Sorting...")
	merged_hits = merge_and_sort_hits( hits_per_file )
	sys.stderr.write("Total sorted hits: {}\n".format(len(merged_hits)))
	write_hits_to_file(merged_hits,args.output_file,args.alignments)


class blastxml_parser():
	def __init__(self,args,sequence_subset):
		self.eval_threshold = args.threshold
		self.include_alignment = args.alignments
		self.max_hits = args.max_hits
		self.sequence_subset = sequence_subset
		self.invert_subset = False

	def __call__(self,blast_filename):
		hits = []
		with open(blast_filename,"r") as blast_fh:
			iterations = NCBIXML.parse(blast_fh)
			formatted_source_file = format_source_file(blast_filename)
			sys.stderr.write("Processing {}\n".format(blast_filename))
			for n_it,iteration in enumerate(iterations):
				iteration_hits = process_iteration(iteration, self.eval_threshold,formatted_source_file, self.include_alignment, self.max_hits) \
									if sequence_should_be_processed(iteration,self.sequence_subset, self.invert_subset) else False
				if iteration_hits:
					hits += iteration_hits
				if n_it % 100000 == 0:
					sys.stderr.write("{} iteration {}\n".format(blast_filename,n_it))
			sys.stderr.write("Finished processing {}, {} hits kept\n".format(blast_filename,len(hits)))
		return hits

def process_iteration(iteration_obj,eval_threshold, source_file, include_alignment,max_hits):
	iteration_hits = []
	for hit_num,hit in enumerate(iteration_obj.alignments):
		hsp = hit.hsps[0] #Only first hsp
		#Store relevant fields in a tuple
		if hsp.expect < eval_threshold:
			visual_alignment = None
			if include_alignment:
				visual_alignment = ( hsp.sbjct_start, hsp.sbjct, hsp.sbjct_end, hsp.match, hsp.query_start,hsp.query,hsp.query_end )
			species = extract_species_name(hit.hit_def)
			pct_id, pct_pos, q_cov = calc_hsp_stats(hsp, iteration_obj.query_length)
			result = (hsp.expect,
					species,
					iteration_obj.query_length,
					q_cov,
					pct_id,
					pct_pos,
					iteration_obj.query,
					source_file,
					hit.hit_def,
					visual_alignment)
			iteration_hits.append(result)

		#Limit number of hits to report
		if hit_num >= (max_hits-1):
			break
	return iteration_hits


def sequence_should_be_processed(iteration_obj, sequence_subset, reverse=False):
	# != in this context corresponds to xor operation!
	return not sequence_subset or ( (iteration_obj.query in sequence_subset) != reverse )

def format_source_file(source_file):
	return os.path.basename(source_file).rsplit(".")[0]

#Ideally annotate from TaxID
def extract_species_name(hit_name):
	species = hit_name
	blastx_species = re.search(r"\[(.+?)\]",species)
	if blastx_species:
		species = blastx_species.group(1)
	else:
		species = hit_name.split(",")[0]
	return species

#Calculates percent identity, percent positives and query coverage
def calc_hsp_stats( hsp , qseq_len):
	pct_id = 100.0 * hsp.identities / hsp.align_length
	pct_pos = 100.0 * hsp.positives / hsp.align_length
	q_cov = 100.0 * (hsp.query_end - hsp.query_start) / qseq_len
	return pct_id, pct_pos, q_cov

def merge_and_sort_hits(hits_per_file):
	combined_hits = [x for x in itertools.chain(*hits_per_file)]
	combined_hits.sort(key=lambda item:item[0],reverse=False)
	return combined_hits

def write_hits_to_file(hits,out_fh,write_alignments=False):
	output_columns = ["e-value", "species", "seq_len", "q_cov","pct_id","pct_pos","ctg_id","source","hit_name"]
	out_fh.write("\t".join(output_columns)+"\n")
	for res in hits:
		out_fh.write( "\t".join([str(x) for x in res[:9]])+"\n")
		if write_alignments:
			#print alignments as well
			alignment = res[9]
			out_fh.write( "\n"+formatAlignment(alignment)+"\n\n" )

# Alignment receives a tuple with the following fields:
# hsp.sbjct_start, hsp.sbjct, hsp.sbjct_end, hsp.match, hsp.query_start,hsp.query,hsp.query_end
def formatAlignment(alignment,space_padding="   "):
	pos_padding = max(len(str(alignment[0])), len(str(alignment[4])) )

	ref_line = space_padding + ("   Ref: {: <"+ str(pos_padding) + "}{}").format( alignment[0] , space_padding )
	qry_line = space_padding + ("   Qry: {: <"+ str(pos_padding) + "}{}").format( alignment[4] , space_padding )
	midline = " " * len(ref_line)

	#Add sequences , padding and end coordinates
	ref_line += alignment[1] + space_padding + str(alignment[2])
	midline  += alignment[3]
	qry_line += alignment[5] + space_padding + str(alignment[6])

	return ref_line+"\n"+midline+"\n"+qry_line

#******************************************************************************
#			Console argument processing
#******************************************************************************

#Assumes files and folder have been validated
def get_files_to_process(args):
	files_to_process = []
	for current_file in filter(os.path.isfile, args.input):
		files_to_process.append(current_file)
	for folder in filter(os.path.isdir, args.input):
		for root_folder,dirs,files in os.walk(folder):
				files_to_process += [os.path.join(root_folder,f) for f in files if f.endswith(".xml") and "blast" in f ]
	#Filters
	if args.filter:
		files_to_process = filter(lambda x:args.filter in x, files_to_process)
	if args.filter_out:
		files_to_process = filter(lambda x: args.filter_out not in x, files_to_process)

	return files_to_process


def load_sequence_subset(filename):
	seq_subset = []
	with open(filename) as seqlist_fh:
		for line in seqlist_fh:
			seq_subset.append(line.rstrip("\n"))
	return frozenset(seq_subset)


def create_argparser():
	parser = argparse.ArgumentParser(description="The program extracts the highest scoring hits for each query \
		sequence and outputs them ordered by e-value. Requires filenames to finish in _tool_blast[nxp]_database.xml")

	parser.add_argument("input",nargs="+",help="Folder in which to look for blast xml files")
	parser.add_argument("-p","--num_proc", type=int, default=1, help="Number of processes to parse the files concurrently. There must be enough memory to keep both in memory")
	parser.add_argument("-o","--output-file", type=argparse.FileType('w'), default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")

	parser.add_argument("--max-hits",type=int,default="10",help="Max number of hits to report for each query sequence")
	parser.add_argument("-a","--alignments",action='store_true',help="Include a visualization for each alignment on the file")
	parser.add_argument("-t","--threshold",type=float,default=1,help="Exclude hits with e-value lower than threshold. Default is '1' (all hits)")

	parser.add_argument("-f","--filter",type=str,help="Only process files that contain the filter string")
	parser.add_argument("-F","--filter-out",type=str,help="Only process files that DO NOT contain the filter string")

	parser.add_argument("-s","--subset-file",type=str,default=None,help="Process only hits from sequences in specified txt file")

	return parser.parse_args()


if __name__ == '__main__':
	main()

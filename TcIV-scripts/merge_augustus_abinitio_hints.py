#!/usr/bin/env python
"""
Script that merges gene annotation from Augustus ab initio and hints modules into a single annotation.
Information from the hints module is priorized and the annotation is completed or corrected with 
ab initio annotations that are not redudant to information from the hints module.

Assumes both files are ordered by contig name.

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
import logging

from collections import OrderedDict,deque
import re

#*************CONSTANTS***********
#GFF column indexes
CTG_ID = 0
ENTRY_TYPE = 2 #gene, transcript, start_codon, stop_codon, CDS
START_POS = 3
STOP_POS = 4
ATTRIBS = 8

#Loaded gff genes
GEN_INTERVAL = 0
GFF_LINES = 1

#****************Begin of Main ***************
def main(args):
	#1) Load data from Ab Initio annotation and EST hints into memory
	abinitio_data = loadgff(args.ab_initio)
	hints_data = loadgff(args.hints)

	#2) Merge both annotations, giving priority to evidence from hints
	current_abinitio = None
	current_hints = None
	abinit_ctg, abinit_genes  = abinitio_data.popitem(False)
	hint_ctg, hint_genes = hints_data.popitem(False)
	#For each contig
	while abinitio_data and hints_data:
		if abinit_ctg == hint_ctg:
			#2.1) If there is both ab initio and hints annotation for the current contig
			#	keep non redundant annotations from both, priorizing hints over ab initio
			mergeGenes(abinit_ctg, abinit_genes,hint_genes, args.output_file)
			abinit_ctg, abinit_genes  = abinitio_data.popitem(False)
			hint_ctg, hint_genes = hints_data.popitem(False)
		
		elif abinit_ctg < hint_ctg:
			#2.2) If there is only ab initio annotation for the current contig
			#	Write all ab initio annotation of that contig
			writeWholeContig(abinit_genes, args.output_file,is_abinitio=True)
			abinit_ctg, abinit_genes  = abinitio_data.popitem(False)
		
		else:
			#2.3) If there is only hints annotation for the current contig
			#	Write all hints annotation of that contig
			writeWholeContig(hint_genes, args.output_file,is_abinitio= False)
			hint_ctg, hint_genes = hints_data.popitem(False)
	
	#3)Verify at least one of the annotation lists has been completely processed		
	assert (len(abinitio_data) == 0) or (len(hints_data) == 0 )
	
	#4) Write remaining annotation 
	while abinitio_data:
		writeWholeContig(abinit_genes, args.output_file,True)
		abinit_ctg, abinit_genes  = abinitio_data.popitem(False)
	
	while hints_data:
		writeWholeContig(hint_genes, args.output_file,False)
		hint_ctg, hint_genes = hints_data.popitem(False)


#*****************End of Main**********************	

#****************Functions**********************
#Extracts the ID attribute from a line, of the format ID=<id>
def extractID(line):
	id_match = re.search(r"ID=(.+?)(;|\n|$)",line)
	gene_id = id_match.group(1) if id_match else None
	return gene_id

#Extracts the Parent attribut from a line, of the format Parent=<parent>
def extractParent(line):
	parent_id_match = re.search(r"Parent=(.+?)(\..+?)?(;|\n|$)",line)
	gene_id = parent_id_match.group(1) if parent_id_match else None
	return gene_id

#Checks if the interval from <contained> is `fully contained` by <container>
#e.g  
#fullycontains( (1,4) , (2,3) ) --> True
def fullyContains(container, contained):
	return (container[0] <= contained[0]) and (contained[1] <= container[1])

def no_overlap(coord1,coord2):
	return coord1[1] < coord2[0] or coord2[1] < coord1[0]

def pct_overlap(main_annot,secondary_annot):
	ovlp_pct = 0
	#if there is some overlap
	if not(main_annot[1] < secondary_annot[0] or secondary_annot[1] < main_annot[0] ):
		ovlp_size = min(main_annot[1],secondary_annot[1]) - max(main_annot[0],secondary_annot[0])
		#Overlap pct defined as the % of the main annotation length covered by the secondary annotation
		ovlp_pct = (ovlp_size*100)/(main_annot[1]-main_annot[0])      
	return ovlp_pct

# Defines the merging rules for both files
# 
# Merging rules:
#
#If there is no overlap:
#	Write the "gene" with the lower coordinates
#If there is one fully contains the other:
#	if the "hints" fully contains the "ab initio"
#		Write hints and discard ab initio
#	if "ab initio" fully contains the "hints" interval
#		Write "ab initio" only if hints covers less than 65% of the length of "ab initio"
#		Write "hints" otherwise
#If there is a partial overlap:
#	if the overlap covers less than 10% of the "hints" length
#		Consider them different, and write the one with smaller coordinates
#	else:
#		Discard the ab initio as "redundant" to the hints annotation		 
def mergeGenes(contig,abinit_genes, hint_genes, merged_fh):
	abinit_coord, abinit_gff = abinit_genes.popleft()
	hint_coord, hint_gff = hint_genes.popleft()

	while abinit_genes and hint_genes:
		if no_overlap(hint_coord, abinit_coord):
			if hint_coord[0] < abinit_coord[0]:
				merged_fh.write(hint_gff)
				hint_coord, hint_gff = hint_genes.popleft()
			else:
				merged_fh.write(abInitioDataDecorator(abinit_gff))
				abinit_coord, abinit_gff = abinit_genes.popleft()

		else: #If there is some overlap
			if fullyContains(hint_coord,abinit_coord): #If hint interval contains the ab initio
				#Keep hints, discard ab initio
				merged_fh.write(hint_gff) #Do I want to write this out? 
				abinit_coord, abinit_gff = abinit_genes.popleft()
				hint_coord, hint_gff = hint_genes.popleft()
			
			elif fullyContains(abinit_coord,hint_coord): #If the ab initio interval contains completely the hints interval
				#If the hits interval spans less the 65% of the ab initio interval - Use ab initio
				if pct_overlap(abinit_coord,hint_coord ) < 65:
					merged_fh.write(abInitioDataDecorator(abinit_gff) )
				else: #If it is covered more than 65% , then keep the hits 
					merged_fh.write(hint_gff)
				abinit_coord, abinit_gff = abinit_genes.popleft()
				hint_coord, hint_gff = hint_genes.popleft()
				
			else: #Overlap is only partial
				#If overlap spans less than 10% of the length of hint, consider them different
				if pct_overlap(hint_coord, abinit_coord) <= 10:
					if hint_coord[0] < abinit_coord[0]:
						merged_fh.write(hint_gff)
						hint_coord, hint_gff = hint_genes.popleft()
					else:
						merged_fh.write(abInitioDataDecorator(abinit_gff))
						abinit_coord, abinit_gff = abinit_genes.popleft()
				else: #If overlap covers more than 10% of hint, discard ab_init
					logging.info("Discarded ab-initio, for more than 10% overlap:"+contig+" "+str(abinit_coord)+"\n")
					abinit_coord, abinit_gff = abinit_genes.popleft()            

#Adds the "src=abinit" attribute at the end of the gff line when the annotation came from ab initio
def abInitioDataDecorator( ab_init_data):
	return re.sub(r"\t(.*)\n",r"\t\g<1>;src=abinit\n",ab_init_data)

#****************I/O functions*******************
def loadgff(gff_filehandle):
	data = OrderedDict()

	current_ctg = ""
	current_gene_id = ""
	current_gene_info = ""
	current_interval = None
	capture_prot_seq = False
	for line in gff_filehandle:
		if line and line != "\n" and line[0] != "#" : #ignore comment lines
			gff_fields = line.rstrip("\n").split("\t")
			
			if gff_fields[ENTRY_TYPE] == "gene":
				if current_gene_id != "":
					#Initialize list for each contig
					if current_ctg not in data:
						data[current_ctg] = deque()
					data[current_ctg].append( (current_interval, current_gene_info)   )
				
				current_ctg = gff_fields[CTG_ID]
				current_gene_id = extractID(gff_fields[ATTRIBS])
				current_interval = ( int(gff_fields[START_POS]) , int(gff_fields[STOP_POS])  )
				current_gene_info = line
			else: #Accumulate under gene
				parent_id = extractParent(gff_fields[ATTRIBS])
				assert( current_gene_id == parent_id)
				current_gene_info += line
		elif line.startswith("# protein sequence"):
			capture_prot_seq = True
			current_gene_info += line
		elif capture_prot_seq:
			current_gene_info += line
			capture_prot_seq = not line.endswith("]\n")
		elif line.rstrip("\n") == "# command line:":
			break

	return data

def writeWholeContig(genes, merged_fh,is_abinitio):
	while genes:
		_,data = genes.pop()
		data = data if not is_abinitio else abInitioDataDecorator(data)
		merged_fh.write(data)


def validate_args(args):
	return True

if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="Merges gene annotation from Augustus ab initio and hints \
		module into a single annotation, priorizing information from the hints module when annotation from both modules are redundant. Assumes both files are lexicographically ordered increasingly by contig name and that that for each contig coordinates are also increasingly ordered")
	parser.add_argument("ab_initio",help="gff file from Augustus ab initio gene prediction", type=file)
	parser.add_argument("hints",help="gff file from Augustus hints based gene prediction", type=file)

	parser.add_argument("-o","--output-file", type=argparse.FileType("w"), default=sys.stdout, help="Name of the output file" )

	parser.add_argument("-l","--log-file", default=None, help="Name of the log file")
	#parser.add_argument("-v","--verbose", action="store_true", help="Print extra information about the execution")
			
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
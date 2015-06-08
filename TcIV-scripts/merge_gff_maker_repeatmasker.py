#!/usr/bin/env python
"""
Script that merges gff files from Maker and Repeat Masker

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

# Time of script execution
import time
import logging

from collections import defaultdict, deque, OrderedDict
import heapq

# *************CONSTANTS***********
# GFF2/3 column indexes
SEQ_ID = 0
ENTRY_TYPE = 2 #gene, transcript, start_codon, stop_codon, CDS
START_POS = 3
STOP_POS = 4
#Only GFF3
GFF3_ATTRIBS = 8


# ****************Begin of Main ***************
def main(argv):
    # 1) Load data from Repeat Masker
    repeatmasker_gff = load_gff(argv.repeatmasker_gff,repeatmasker_gff_iterator)

    # 2) Load data from Maker
    maker_gff = load_gff(argv.maker_gff,maker_gff_iterator)

    rm_ctg, rm_genes = repeatmasker_gff.popleft()
    maker_ctg, maker_genes = maker_gff.popleft()

#    while repeatmasker_gff and maker_gff:
    while len(repeatmasker_gff) and len(maker_gff):
        if rm_ctg == maker_ctg:
            # 2.1) If there is both ab initio and hints annotation for the current contig
            # keep non redundant annotations from both, priorizing hints over ab initio
            merge_features(maker_ctg,maker_genes, rm_genes, argv.output_file)
            rm_ctg, rm_genes  = repeatmasker_gff.popleft()
            maker_ctg, maker_genes = maker_gff.popleft()

        elif rm_ctg < maker_ctg:
            # 2.2) If there is only repeat masker annotation for the current contig
            write_all_genes_in_contig(rm_genes, argv.output_file,True)
            rm_ctg, rm_genes  = repeatmasker_gff.popleft()

        else:
            # 2.3) If there is only maker annotation for the current contig
            write_all_genes_in_contig(maker_genes, argv.output_file,False)
            maker_ctg, maker_genes = maker_gff.popleft()

    # 3)Verify at least one of the annotation lists has been completely processed
    assert (len(repeatmasker_gff) == 0) or (len(maker_gff) == 0 )

    # 4) Write remaining annotation
    while repeatmasker_gff:
        write_all_genes_in_contig(rm_genes, argv.output_file,True)
        rm_ctg, rm_genes  = repeatmasker_gff.popitem(False)

    while maker_gff:
        write_all_genes_in_contig(maker_genes, argv.output_file,False)
        maker_ctg, maker_genes = maker_gff.popitem(False)


# End of Main *****************************************
def merge_features(seq_id,main_features,secondary_features,out_fh):
    main_pos, main_feat = main_features.popleft()
    secondary_pos, secondary_feat = secondary_features.popleft()

    while main_features and secondary_features:
        if fully_contains(main_pos,secondary_pos) or fully_contains(secondary_pos,main_pos):
            out_fh.write(main_feat)
            main_pos, main_feat = main_features.popleft()
            secondary_pos, secondary_feat = secondary_features.popleft()
            if fully_contains(secondary_pos,main_pos):
                sys.stderr.write("WARNING: Secondary covers more than primary in {}:{}\n".format(seq_id,main_pos))
        else:
            if not no_overlap(main_pos,secondary_pos):
                sys.stderr.write("WARNING: {}% of overlap in {}:{}\n".format(pct_overlap(main_pos,secondary_pos),seq_id,main_pos))
            if main_pos < secondary_pos:
                out_fh.write(main_feat)
                main_pos, main_feat = main_features.popleft()
            else:
                out_fh.write(secondary_feat)
                secondary_pos, secondary_feat = secondary_features.popleft()

    while main_features:
        out_fh.write(main_feat)
        _, main_feat = main_features.popleft()

    while secondary_features:
        out_fh.write(main_feat)
        _, secondary_feat = secondary_features.popleft()



def write_all_genes_in_contig(genes, out_fh,is_secondary):
    while genes:
        _, data = genes.pop()
        data = data if not is_secondary else secondaryDataDecorator(data)
        out_fh.write(data)


def secondaryDataDecorator(data):
    return data

# *****************End of Main**********************

# Checks if the interval from <contained> is `fully contained` by <container>
#   e.g fullycontains( (1,4) , (2,3) ) --> True
def fully_contains(container, contained):
    return (container[0] <= contained[0]) and (contained[1] <= container[1])


def no_overlap(coord1,coord2):
    return coord1[1] <= coord2[0] or coord2[1] <= coord1[0]

def pct_overlap(main_annot,secondary_annot):
    ovlp_pct = 0
    # if there is some overlap
    if not(main_annot[1] < secondary_annot[0] or secondary_annot[1] < main_annot[0] ):
        ovlp_size = min(main_annot[1],secondary_annot[1]) - max(main_annot[0],secondary_annot[0])
        # Overlap pct defined as the % of the main annotation length covered by the secondary annotation
        ovlp_pct = (ovlp_size*100.0)/(main_annot[1]-main_annot[0])
    return ovlp_pct


# ****************I/O functions*******************
def load_gff(gff_filehandle,gff_file_iterator):
    data = defaultdict(deque)
    for ctg, pos, gene in gff_file_iterator(gff_filehandle):
        # Push data into deque
        data[ctg].append((pos, gene))
    ordered_data = deque((ctg, data[ctg]) for ctg in sorted(data))
    return ordered_data


def maker_gff_iterator(maker_gff_fh):
    current_gene_info = ""
    current_ctg = ""
    for line in maker_gff_fh:
        if line and line != "\n" and line[0] != "#":  # ignore comment lines
            gff_fields = line.rstrip("\n").split("\t")

            if gff_fields[ENTRY_TYPE] == "gene":
                if current_gene_info != "":
                    yield (current_ctg, current_pos, current_gene_info)

                current_ctg = gff_fields[SEQ_ID]
                current_pos = (int(gff_fields[START_POS]), int(gff_fields[STOP_POS]))
                current_gene_info = line
            else:  # Accumulate under gene
                current_gene_info += line

    yield (current_ctg,current_pos,current_gene_info)


def repeatmasker_gff_iterator(repeatmasker_gff_fh):
    for line in repeatmasker_gff_fh:
        if line.startswith("#"):
            continue
        gff_fields = line.rstrip("\n").split("\t")
        current_ctg = gff_fields[SEQ_ID]
        current_pos = (int(gff_fields[START_POS]), int(gff_fields[STOP_POS]))
        current_gene_info = format_repeatmasker_line(line)
        yield (current_ctg,current_pos,current_gene_info)


def format_repeatmasker_line( rm_line ):
    return rm_line.replace("similarity","repeat",1).replace("Target ","id=")


def validate_args(argv):
    return True

if __name__ == '__main__':
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Completes annotation from MAKER with detected repeats with RepeatMasker")
    parser.add_argument("maker_gff",metavar="maker-gff",help="gff file from Augustus ab initio gene prediction", type=file)
    parser.add_argument("repeatmasker_gff",metavar="repeatmasker-gff",help="gff file from Augustus hints based gene prediction", type=file)

    parser.add_argument("-o","--output-file", type=argparse.FileType("w"), default=sys.stdout, help="Name of the output file" )
    parser.add_argument("-l","--log-file", default=None, help="Name of the log file")
    #parser.add_argument("-v","--verbose", action="store_true", help="Print extra information about the execution")

    args = parser.parse_args()

    if validate_args(args):
        # Initialize log
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
#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

import collections

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from matplotlib import collections as mc



#Define Coords 'class' 
Coords = collections.namedtuple('Coords', ['start', 'end'])


def Coords_to_segment(coord,y_pos):
    return [(coord.start,y_pos),(coord.end,y_pos)]


def read_blast_cols(blast_cols_file):
    with open(blast_cols_file) as fh:
        blast_cols = next(fh).rstrip(" \n").split(" ")
    return blast_cols


'''
Assuming 
* reference sequence is the subject sequence
* the reference coordinates goes from 0 to length of sequence

Calculate the coordinates for the aligned query sequences, so that
it would reflect the actual alignment, where the overlapping segments
actually overlap 
'''
def calculate_query_coords_wrt_subject(df):
    #Calculate lines of query hits, shifted with respect to their alignment to subject sequence
    full_queries = []
    overlapping_regions = []
    labels = []

    for _, row in df.iterrows():
        #Positive strand
        if row["sstart"] < row["send"]:
            shifted_start = row["sstart"] - row["qstart"]
        #Negative strand
        else:
            shifted_start = row["sstart"] - row["qend"]

        #Full length of the query sequence, with coords relative to the alignment to the subject
        new_full_query = Coords(shifted_start, shifted_start+row["qlen"])
        #The overlap is always based on the subject coordinates
        new_overlap = Coords(row["sstart"], row["send"])

        full_queries.append(new_full_query)
        overlapping_regions.append(new_overlap)
        labels.append(row["qseqid"])

    return (full_queries, overlapping_regions, labels)

def plot_hits(df, sseqid, out_file): 
    ref_seq_length = df.slen.iloc[0]
    #Calculate reference line (subject sequence)
    reference_seq = Coords_to_segment(Coords(0, ref_seq_length-1), y_pos=0)

    query_seqs, ovlp_regions, labels = calculate_query_coords_wrt_subject(df)

    query_seqs_lines = [Coords_to_segment(seq, (hit_idx+1)*0.1) for hit_idx, seq in enumerate(query_seqs)]
    ovlp_region_lines = [Coords_to_segment(seq, (hit_idx+1)*0.1) for hit_idx, seq in enumerate(ovlp_regions)]

    pad_left_down = [(-ref_seq_length/4, -len(query_seqs)*0.2), (0, -len(query_seqs)*0.2)]
    pad_right_up = [(ref_seq_length, (1+len(query_seqs))*0.1),
                    (5*ref_seq_length/4, (1+len(query_seqs))*0.1)]

    #Matplotlib plot
    mpl.rcParams['figure.figsize'] = (10.0, 8.0)
    fig,ax = plt.subplots()
    lc1 = mc.LineCollection([reference_seq], linewidths=7, colors=["blue"])
    lc2 = mc.LineCollection(query_seqs_lines, linewidths=3, colors=["black"])
    lc3 = mc.LineCollection(ovlp_region_lines, linewidths=3, colors=["green"])
    pad = mc.LineCollection([pad_left_down, pad_right_up], linewidths=3, colors=["white"])

    ax.add_collection(lc1)
    ax.add_collection(lc2)
    ax.add_collection(lc3)
    ax.add_collection(pad)

    #Add labels
    for hit_idx, (lbl, query_seqs) in enumerate(zip(labels, query_seqs)):
        ax.text(query_seqs.start, ((hit_idx+1)*0.1)+0.01, lbl)

    #Set title
    ax.autoscale()
    ax.set_title(sseqid)

    fig.savefig(out_file)
    plt.close()

'''
Main function
'''
def main():
    #Read blast column names
    blast_cols = read_blast_cols(sys.argv[2])

    #Read tsv file
    blast_hits = pd.read_csv(sys.argv[1], header=None, names=blast_cols, sep="\t")

    for sseqid, hits_df in blast_hits.groupby("sseqid"):
        plot_hits(hits_df, sseqid, "{}.png".format(sseqid))

if __name__ == '__main__':
    main()



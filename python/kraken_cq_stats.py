#!/usr/bin/env python
"""
Calculates stats

lolol program
@author: Mauricio Barrientos-Somarribas
"""
from functools import reduce
import operator
import sys


def main():
    """Main method."""
    if len(sys.argv) != 3:
        print("{} <kraken_out> <kraken_db>".format(sys.argv[0]),
              file=sys.stderr)
        exit(1)

    kraken_out = sys.argv[1]
    kraken_db = sys.argv[2]

    tax_tree = load_taxonomy("{}/taxonomy/nodes.dmp".format(kraken_db))

    read_count = 0
    classified_count = 0
    with open(kraken_out) as kraken_out_fh:
        for line in kraken_out_fh:
            read_status, kmer_stats = parse_line(line, tax_tree)
            read_count += 1
            if read_status == "C":
                classified_count += 1
                print(format_line_stats(kmer_stats))

    print("{}\t{}\t{}".format(kraken_out, read_count,
                              classified_count), file=sys.stderr)


def load_taxonomy(nodes_file):
    """Load the nodes.dmp file into a dictionary."""
    tax_tree = {}
    with open(nodes_file, "r") as nodes_fh:
        for line in nodes_fh:
            node_id, parent_id, _ = line.rstrip("\n").split("\t|\t", 2)
            tax_tree[int(node_id)] = int(parent_id)

    tax_tree[1] = 0
    return tax_tree


def parse_line(line, tax_tree):
    """
    Extract relevant information from an entry from kraken output.

    Fields:
    C/U: Classified unclassified
    seqid
    taxid
    seq length
    kmer mappings
    """
    status, _, assigned_taxid, _, kmer_str = line.rstrip("\n").split("\t")

    # Calculate stats only for classified reads
    kmer_stats = None
    if status == "C":
        kmer_stats = calculate_kmer_stats(int(assigned_taxid),
                                          kmer_str, tax_tree)

    return (status, kmer_stats)


def get_tax_branch(taxid, tax_tree):
    """
    Return a list of all parents of a given taxid.

    The list begins with the taxid and ends with the root of the tree.
    """
    branch = []
    while taxid != 0:
        branch.append(taxid)
        # Gets the parent of current tax id
        taxid = tax_tree[taxid]
    return branch


def calculate_kmer_stats(assigned_taxid, kmer_str, tax_tree):
    """Calculate C/Q for assigned taxid and parents and other stats."""
    hits = kmer_str.split(" ")
    taxid_counts = {}
    # 1. Count total of kmers assigned to a given taxid
    # Note, keeps the taxid as a string because of 'A' ambiguous kmers
    for hit in hits:
        taxid, kmer_count = hit.split(":")
        if taxid not in taxid_counts:
            taxid_counts[taxid] = 0
        taxid_counts[taxid] += int(kmer_count)


    stats_to_report = {
        "total_kmers": reduce(operator.add,
                              [taxid_counts[x] for x in taxid_counts]),
        "unclassified_kmers": taxid_counts['0'] if '0' in taxid_counts else 0,
        "ambiguous_kmers": taxid_counts["A"] if 'A' in taxid_counts else 0
    }

    # 2. Calculate total kmer support for assigned taxid and its parent nodes
    #    Count all kmers 'consistent' with the given taxid
    #    i.e the count should be kmers assigned to itself or any parent node

    # Convert to string since `taxid_counts` stores the taxids as 'strings'
    tax_branch = [str(x) for x in get_tax_branch(assigned_taxid, tax_tree)]
    cumulative_count = 0
    for node in tax_branch:
        if node in taxid_counts:
            taxid_counts[node] += cumulative_count
            cumulative_count = taxid_counts[node]

    cum_counts_to_report = []
    for node in tax_branch:
        if node in taxid_counts:
            cum_counts_to_report.append(taxid_counts[node])
        else:
            cum_counts_to_report.append("*")

    # Add to stats report
    stats_to_report["kmer_support"] = cum_counts_to_report

    return stats_to_report


def format_line_stats(kmer_stats):
    """Create a parsable string for the calculated kmer stats."""
    fields = ("total_kmers",
              "unclassified_kmers",
              "ambiguous_kmers",
              "kmer_support")

    formatted_str = "\t".join(str(kmer_stats[x]) for x in fields[:3])
    formatted_str += "\t" + "\x1F".join(str(x) for x in kmer_stats["kmer_support"])
    return formatted_str


if __name__ == '__main__':
    main()

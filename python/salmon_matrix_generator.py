#!/usr/bin/env python
"""
Script that takes an arbitrary number of Salmon output files and generates two matrixes:
(Sample x TPM) and (Sample x NumReads) with the abundance estimates for each gene for all samples


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

__author__ = 'maubarsom'

import sys
import argparse
import pandas as pd

def main(argv):
    sample_names = [ extract_samplename(f) for f in argv.salmon_files ]
    print "Loaded {}".format(" ".join(sample_names))
    salmon_files_list = [pd.read_csv(f,sep='\t',skiprows=12,header=None,names=["Name","Length","TPM","NumReads"])
                    for f in argv.salmon_files]

    tpm_matrix = pd.concat(map(lambda x: x["TPM"], salmon_files_list), axis=1)
    numreads_matrix = pd.concat(map(lambda x: x["NumReads"], salmon_files_list), axis=1)

    #Add gene names column to matrixes
    tpm_matrix = pd.concat([ salmon_files_list[0]["Name"],tpm_matrix ],axis=1)
    numreads_matrix = pd.concat([ salmon_files_list[0]["Name"],numreads_matrix ],axis=1)

    #Add sample names to columns
    tpm_matrix.columns = ["Name"] + sample_names
    numreads_matrix.columns = ["Name"] + sample_names

    tpm_matrix.to_csv("{}_tpm.csv".format(argv.output_prefix),index=False)
    numreads_matrix.to_csv("{}_numreads.csv".format(argv.output_prefix),index=False)

def extract_samplename(salmon_file):
    samplename = None
    with open(salmon_file) as fh:
        for line_number, line in enumerate(fh):
            if line_number == 9: # 10th line
                samplename = line.split("] => {")[1].lstrip(" ").rstrip("\n }")
                break
    return samplename

if __name__ == '__main__':
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Script that creates (TPM/NumReads x Sample) matrixes from Salmon files")
    parser.add_argument("salmon_files",help="Input file",nargs="+")
    parser.add_argument("-o","--output-prefix", default="out", help="Prefix for the output files" )
    args = parser.parse_args()
    main(args)

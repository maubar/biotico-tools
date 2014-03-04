#!/usr/bin/env python
"""
Docstring!


Author: Mauricio Barrientos
Email:  mauricio.barrientos@ki.se

License notice
"""


import sys
import argparse
import os.path

#Time of script execution
from datetime import datetime

import re
import math

import itertools
from collections import * 

#Data Analysis libs
import numpy as np
import matplotlib as plt
import pandas as pd
import scipy.stats as stats


#****************Begin of Main ***************
def main(args):
	pass


#*****************End of Main**********************	
	
#@TODO: Write generic lazy file reader stub
def lazyFileReader(filename):
	pass


def validate_args(args):
	if not os.path.isfile(args.input):
		sys.stderr.write("Input file does not exist!\n")
		sys.exit(1)

	return True


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="WRITE DESCRIPTION HERE",epilog= "TEXT AFTER THE HELP")
	parser.add_argument("input",help="Input file",nargs="?", type=file, default=sys.stdin)
	
	parser.add_argument("-o","--output-file", default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-o","--output-prefix", default=None, help="Prefix of the output file(s)" )
	parser.add_argument("-o","--output-folder", default="./", help="Folder where to output the results" )

	parser.add_argument("-l","--log-file", help="Name of the log file",type=argparse.FileType("w"), default=sys.stderr)
	groupV.add_argument("-v","--verbose", help="Print extra information about the execution",action="store_true")
			
	args = parser.parse_args()
	
	if validate_args(args):
		main( args )
	else:
		sys.stderr.write("Invalid arguments. Exiting script\n")
		sys.exit(1)
#!/usr/bin/env python
"""
Docstring!


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
	return True


if __name__ == '__main__':
	#Process command line arguments
	parser = argparse.ArgumentParser(description="WRITE DESCRIPTION HERE",epilog= "TEXT AFTER THE HELP")
	parser.add_argument("input",help="Input file",nargs="?", type=file, default=sys.stdin)
	
	parser.add_argument("-o","--output-file", type=file, default=sys.stdout, help="Name of the output file" )
	parser.add_argument("-o","--output-prefix", default=None, help="Prefix of the output file(s)" )
	parser.add_argument("-o","--output-folder", default="./", help="Folder where to output the results" )

	parser.add_argument("-l","--log-file",type=argparse.FileType("w"), default=sys.stderr, help="Name of the log file")
	parser.add_argument("-v","--verbose", action="store_true", help="Print extra information about the execution")
			
	args = parser.parse_args()
	
	if validate_args(args):
		main( args )
	else:
		sys.stderr.write("Invalid arguments. Exiting script\n")
		sys.exit(1)
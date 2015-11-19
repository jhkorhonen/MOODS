#!/usr/bin/python

import MOODS.scan
import MOODS.tools
import MOODS.parsers

import math
import sys
import time
import os

# actual program starts here

if len(sys.argv) < 3:
	print "usage: python ex-mixing-data.py [sequence file] [matrix directory]"
	sys.exit(1)

# read an sequence
# in this example we do not do any fancy input parsing,
# so the file should have only characters ACGT
# (so no line breaks for instance)
file_name = sys.argv[1]
with open(file_name, "r") as file_handle:
    seq = file_handle.read()

# read the matrix files
# we'll mix both 0-order JASPAR models and first-order models together
matrix_directory = sys.argv[2]
pfms = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']
adms = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.adm']

# use flat bg for the models
bg = MOODS.tools.flat_bg(4)

# log-odds
lo_pfms = [MOODS.parsers.pfm_log_odds(matrix_directory + filename, bg, 1) for filename in pfms]
lo_adms = [MOODS.parsers.adm_log_odds(matrix_directory + filename, bg, 0.0001) for filename in adms]
matrices =  lo_pfms + lo_adms

matrix_names = pfms + adms
# scanning
# we'll want about 10000 best hits for each matrix, this tries to pick an appropriate threshold by magic
results = MOODS.scan.scan_best_hits_dna(seq, matrices, 10000)

for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    print matrix_name + "|" + str(len(result))

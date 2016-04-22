#!/usr/bin/env python

import MOODS.scan
import MOODS.tools
import MOODS.parsers

import math
import sys
import time
import os

# This example shows how to use the scan_best_hits function to obtain a desired amount
# of highest-scoring hits.

if len(sys.argv) < 3:
	print "usage: python ex-best-hits.py [sequence file] [matrix directory]"
	sys.exit(1)

# read an sequence
file_name = sys.argv[1]
with open(file_name, "r") as file_handle:
    first = file_handle.next().strip()
    if first[0] == '>':
        first = ''
    seq = "".join([first] + [line.strip() for line in file_handle])

# read the matrix files
# we'll mix both 0-order JASPAR models and first-order models together
matrix_directory = sys.argv[2]
pfms = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']
adms = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.adm']

# use flat bg for the models
bg = MOODS.tools.flat_bg(4)
pseudocount = 0.001
lo_pfms = [MOODS.parsers.pfm_to_log_odds(matrix_directory + filename, bg, pseudocount) for filename in pfms]
lo_adms = [MOODS.parsers.adm_to_log_odds(matrix_directory + filename, bg, pseudocount) for filename in adms]
matrices =  lo_pfms + lo_adms

matrix_names = pfms + adms

# scanning

# we'll want about 1000 best hits for each matrix, this tries to pick an appropriate threshold by magic

results = MOODS.scan.scan_best_hits_dna(seq, matrices, 1000)

# Specifically, MOODS tries to guess a good threshold to produce 10000-30000 hits; if this does not succeed, MOODS uses
# binary search to try and refine the threshold. By default we'll do 10 iterations of this, after which we
# return *some* amount of hits that is less than 100000 (may be 0 in some cases)
# 
# Extra parameters can be used to adjust these values:
#
# results = MOODS.scan.scan_best_hits_dna(seq, matrices, target, iterations = 10, MULT = 3, UPPER_MULT = 10, window_size = 7)
# 
# where
#   target       is the number of hits we want
#   iterations   is the number of iterations (0 for unlimited, is guaranteed to finish at some point but may be slow)
#   MULT         MULT*target is the upper bound for the number of hits MOODS *tries* to find (i.e. desired amount of hits is)
#                in the range [target, MULT*target)
#   UPPER_MULT   UPPER_MULT*target is the upper bound for the number of hits MOODS is allowed to give (safeguards against slow
#                performance from too many hits)
#   window_size  window size for scanning, passed to scanner construction (7 is probably fine always)
#
# Currently you can omit some of these parameters to use defaults, but you can only omit the last ones. For example, this does the
# same thing as above but without iteration limit:
#
# results = MOODS.scan.scan_best_hits_dna(seq, matrices, 1000, 0)



for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    print matrix_name + "|" + str(len(result))

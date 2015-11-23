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

# estimate background from the sequence file
bg = MOODS.tools.bg_from_sequence_dna(seq,1)

# log-odds
lo_pfms = [MOODS.parsers.pfm_log_odds(matrix_directory + filename, bg, 1) for filename in pfms]
lo_adms = [MOODS.parsers.adm_log_odds(matrix_directory + filename, bg, 0.0001) for filename in adms]
matrices =  lo_pfms + lo_adms

matrix_names = pfms + adms

# threshold selection
th = 0.00001
thresholds = [MOODS.tools.threshold_from_p(m,bg,th) for m in lo_pfms] + [MOODS.tools.threshold_from_p(m,bg,th,4) for m in lo_adms]

# scanning
results = MOODS.scan.scan_dna(seq, matrices, bg, thresholds, 7)

for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    if matrix_name in pfms:
        l = len(matrix[0])
    else:
        l = len(matrix[0]) + 1
    for r in sorted(result, key=lambda r: r.pos):
        hitseq = seq[r.pos:r.pos+l]
        print matrix_name + '|' + str(r.pos) + '|' + hitseq + '|'  + str(r.score)

print ""
for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    print matrix_name + "|" + str(len(result))

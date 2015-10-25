#!/usr/bin/python

import MOODS.scan
import MOODS.tools

import math
import sys
import time
import os


# helper function, loads JASPAR matrix files (.pfm)
def load_jaspar(filename):
    with open(filename, "r") as file_handle:
        mat = [[float(i) for i in line.split()] for line in file_handle]
        return mat

# helper function, loads .adm files
# returns first and last zero-order terms for log-odds computation
def load_adm(filename):
    with open(filename, "r") as file_handle:
        lines = [[i for i in line.split()] for line in file_handle]
        mat = [[float(i) for i in line[:-2]] for line in lines[:16]]
        zero_terms = [float(line[0]) for line in lines[-4:]]
        return (mat, zero_terms)

# MOODS does not have log-odds computation for higher-order matrices yet,
# so we'll use this quick hack
def log_odds_dna_ho(m, bg):
	mat = m[0]
	zero_terms = m[1]
	los = [[0] * len(row) for row in mat]

	for i in xrange(len(mat[0])):
		for j in xrange(16):
			los[j][i] = math.log(mat[j][i]) - math.log(bg[j & 3])

	for k in xrange(4):
		for j in xrange(4):
			los[(j << 2) ^ k][0] = los[(j << 2) ^ k][0] + math.log(zero_terms[j]) + math.log(bg[j])
	return los



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

mz = [load_jaspar(matrix_directory + filename) for filename in pfms]
mf = [load_adm(matrix_directory + filename) for filename in adms]

# estimate background from the sequence file
bg = MOODS.tools.bg_from_sequence_dna(seq,1)

# log-odds
lo_pfms = [MOODS.tools.log_odds(m, bg, 1) for m in mz]
lo_adms = [log_odds_dna_ho(m, bg) for m in mf]
matrices =  lo_pfms + lo_adms

matrix_names = pfms + adms

# threshold selection
th = 0.000001
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

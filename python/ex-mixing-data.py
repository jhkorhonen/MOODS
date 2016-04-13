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
# bg = MOODS.tools.flat_bg(4)

# log-odds
lo_pfms = [MOODS.parsers.pfm_log_odds(matrix_directory + filename, bg, 1) for filename in pfms]
# lo_adms = [MOODS.parsers.adm_log_odds(matrix_directory + filename, bg, 0.0001) for filename in adms]
lo_adms = []
matrices =  lo_pfms + lo_adms

matrix_names = pfms + adms

# threshold selection
th = 0
thresholds = [MOODS.tools.threshold_from_p(m,bg,th) for m in lo_pfms] + [MOODS.tools.threshold_from_p(m,bg,th,4) for m in lo_adms]


def consensus(mat):
    ds = {0: "A", 1: "C", 2: "G", 3: "T"}
    ls = []
    for i in range(len(mat[0])):
        row = [line[i] for line in mat]
        ls.append(row.index(max(row)))
    return "".join([ds[x] for x in ls])

def near_max(mat):
    ls = []
    ls.append(mat[0][0])
    for i in range(1,len(mat[0])):
        row = [line[i] for line in mat]
        ls.append(max(row))
    return sum(ls)

# thresholds = [MOODS.tools.max_score(m)-0.01 for m in lo_pfms] + [MOODS.tools.max_score(m,4)-0.01 for m in lo_adms]
# thresholds = [near_max(m) for m in lo_pfms]


for (matrix,matrix_name,threshold) in zip(matrices, matrix_names, thresholds):
    print matrix_name, "(", len(matrix[0]) ,"):", threshold, "/", MOODS.tools.max_score(matrix), consensus(matrix)

# print seq.count("CAATTATT")

# scanning
results = MOODS.scan.scan_dna(seq, matrices, bg, thresholds, 7)
#
# for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
#     if matrix_name in pfms:
#         l = len(matrix[0])
#     else:
#         l = len(matrix[0]) + 1
#     for r in sorted(result, key=lambda r: r.pos):
#         hitseq = seq[r.pos:r.pos+l]
#         print matrix_name + '|' + str(r.pos) + '|' + hitseq + '|'  + str(r.score)
#
print ""
for (matrix,matrix_name,result,th) in zip(matrices, matrix_names, results,thresholds):
    print matrix_name + "|" + str(len(result)) + "|" + str(seq.count(consensus(matrix))) + "|" + str(len(MOODS.scan.naive_scan_dna(seq,matrix,MOODS.tools.max_score(matrix))))
    # if matrix_name in ["MA0001.pfm", "MA0010.pfm"]:
    #     for res in result:
    #         print seq[res.pos:res.pos+len(matrix[0])], res.score
    #     print "--"
    #     for res in MOODS.scan.naive_scan_dna(seq,matrix,MOODS.tools.max_score(matrix)):
    #         print seq[res.pos:res.pos+len(matrix[0])], res.score

    # print matrix_name + "|" + str(len(result)) + "|" + str(len(MOODS.scan.naive_scan_dna(seq,matrix,th)))
    

print ""
#
# for (matrix, matrix_name,th) in zip(matrices, matrix_names,thresholds):
#     print matrix_name, th, len(matrix[0]), consensus(matrix)
#     result = MOODS.scan.scan_dna(consensus(matrix), [matrix], bg, [th], 7)
#     result = MOODS.scan.scan_dna(consensus(matrix), [matrix], bg, [MOODS.tools.min_score(matrix)], 7)
#     print matrix_name + "|" + str(result[0][0].score)
#     print " "*len(matrix_name + "|") + str(th)
    

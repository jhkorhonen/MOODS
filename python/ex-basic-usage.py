#!/usr/bin/python

import MOODS.scan
import MOODS.tools

import sys
import time
import os

def load_matrix(filename):
    with open(filename, "r") as file_handle:
        mat = [[float(i) for i in line.split()] for line in file_handle]
        return mat


total = time.clock()
start = time.clock()

if len(sys.argv) < 3:
	print "usage: python ex-basic-usage.py [sequence file] [matrix directory]"
	sys.exit(1)

# read an sequence
# in this example we do not do any fancy input parsing,
# so the file should have only characters ACGT
# (so no line breaks for instance)
file_name = sys.argv[1]
with open(file_name, "r") as file_handle:
    seq = file_handle.read()

end = time.clock()

print "Reading sequence:", end - start

start = time.clock()

# read the matrix files
# should be in JASPAR .pfm format
matrix_directory = sys.argv[2]
matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']

matrices = [load_matrix(matrix_directory + filename) for filename in matrix_names]
# reverse complements
rc_matrices = [MOODS.tools.reverse_complement(m) for m in matrices]

end = time.clock()

print "Reading matrices:", end - start

start = time.clock()

# estimate background from the sequence file
bg = MOODS.tools.bg_from_sequence_dna(seq,1)

end = time.clock()

print "Computing background:", end - start

# log-odds transformation for matrices
log_odds_matrices = [MOODS.tools.log_odds(m, bg, 1) for m in matrices] + [MOODS.tools.log_odds(m, bg, 1) for m in rc_matrices]
# thresholds computed from p-value
thresholds = [MOODS.tools.threshold_from_p(m, bg, 0.00001) for m in log_odds_matrices]

start = time.clock()

# scanning
results = MOODS.scan.scan_dna(seq, log_odds_matrices, bg, thresholds, 7)

# separate reverse complements and the non-reverse complements 
fr = results[:len(matrices)]
rr = results[len(matrices):]

# mix the results together, use + and - to indicate strand
results = [ [(r.pos, r.score, '+') for r in fr[i]] + [(r.pos, r.score, '-') for r in rr[i]] for i in xrange(len(matrices))]


end = time.clock()

print "Scanning:", end - start

print "Total hits:", sum([len(r) for r in results])

for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    l = len(matrix[0])
    for r in sorted(result, key=lambda r: r[0]):
        strand = r[2]
        pos = r[0]
        hitseq = seq[pos:pos+l]
        print matrix_name + '|' + str(pos) + '|' + strand + '|' + hitseq + '|'  + str(r[1])

# for i, m in zip(range(len(matches)), matches):
#     print "Hits for", i, ":", len(m)

print "Total time:", time.clock() - total
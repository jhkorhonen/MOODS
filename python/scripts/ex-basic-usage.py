#!/usr/bin/env python

import MOODS.scan
import MOODS.tools
import MOODS.parsers

import sys
import time
import os

total = time.clock()
start = time.clock()

if len(sys.argv) < 3:
    print "usage: python ex-basic-usage.py [sequence file] [matrix directory]"
    sys.exit(1)


# ---- parameters ----
pseudocount = 0.00001
pvalue = 0.0001
# ---- end parameters ----


# ---- read sequence ----
# read an sequence
# should be simple single-sequence fasta/plain text
file_name = sys.argv[1]
with open(file_name, "r") as file_handle:
    first = file_handle.next().strip()
    if first[0] == '>':
        first = ''
    seq = "".join([first] + [line.strip() for line in file_handle])

end = time.clock()

print "Reading sequence:", end - start
print "Sequence lenght:", len(seq)
# ---- end read sequence ----


# ---- compute bg ----
start = time.clock()

# estimate background from the sequence file
bg = MOODS.tools.bg_from_sequence_dna(seq,1)

end = time.clock()

print "Computing background:", end - start
# ---- end compute bg ----


# ---- process matrices ----
start = time.clock()

# read the matrix files
# should be in JASPAR .pfm format
matrix_directory = sys.argv[2]
matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']

matrices = [MOODS.parsers.pfm_to_log_odds(matrix_directory + filename, bg, pseudocount) for filename in matrix_names]
# reverse complements
matrices = matrices + [MOODS.tools.reverse_complement(m) for m in matrices]

thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]


end = time.clock()

print "Reading matrices:", end - start
# ---- end process matrices ----


# ---- scanning ----
start = time.clock()
results = MOODS.scan.scan_dna(seq, matrices, bg, thresholds, 7)
end = time.clock()

print "Scanning:", end - start

# ---- end scanning ----

# ---- process results ----
# separate reverse complements and the non-reverse complements 
fr = results[:len(matrix_names)]
rr = results[len(matrix_names):]

# mix the results together, use + and - to indicate strand
results = [ [(r.pos, r.score, '+') for r in fr[i]] + [(r.pos, r.score, '-') for r in rr[i]] for i in xrange(len(matrix_names))]

for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    l = len(matrix[0])
    for r in sorted(result, key=lambda r: r[0]):
        strand = r[2]
        pos = r[0]
        hitseq = seq[pos:pos+l]
        print matrix_name + '|' + str(pos) + '|' + strand + '|' + hitseq + '|'  + str(r[1])

print "Total hits:", sum([len(r) for r in results])
for i, name, m in zip(range(len(results)), matrix_names, results):
    print "Hits for", name, ":", len(m)

print "Total time:", time.clock() - total

# ---- end process results ----
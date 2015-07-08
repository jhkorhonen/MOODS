import MOODS.scan
import MOODS.tools

import time
import os

def load_matrix(filename):
    with open(filename, "r") as file_handle:
        mat = [[float(i) for i in line.split()] for line in file_handle]
        return mat


total = time.clock()
start = time.clock()

# with open("../examples/data/sequence/dnaACGT.txt", "r") as file_handle:
with open("acgt_random_60.txt", "r") as file_handle:
    seq = file_handle.read()
    # seq = ''.join([filter(lambda x: x in 'acgtACGT', line) for line in file_handle if line[0] != '>'])

end = time.clock()

print "Reading sequence:", end - start

start = time.clock()

matrix_directory = "../examples/data/matrix/JASPAR_CORE_2008/"
matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']

matrices = [load_matrix(matrix_directory + filename) for filename in matrix_names]

end = time.clock()

print "Reading matrices:", end - start

start = time.clock()

bg = MOODS.tools.bg_from_sequence_dna(seq,1)

end = time.clock()

print "Computing background:", end - start

matrices = [MOODS.tools.log_odds(m, bg, 1) for m in matrices]
thresholds = [MOODS.tools.threshold_from_p(m, bg, 0.0001) for m in matrices]

results = MOODS.scan.scan_dna(seq, matrices, bg, thresholds, 7)

# for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
#     print matrix_name, len(result)

print "Hits:", sum([len(r) for r in results])

# for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
#     l = len(matrix[0])
#     for r in sorted(result, key=lambda r: r.pos):
#         strand = "+"
#         hitseq = seq[r.pos:r.pos+l]
#         print matrix_name + '|' + str(r.pos) + '|' + strand + '|' + hitseq + '|'  + str(r.score)

# for i, m in zip(range(len(matches)), matches):
#     print "Hits for", i, ":", len(m)

print "Total time for everything:", time.clock() - total
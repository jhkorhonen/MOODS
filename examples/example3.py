import MOODS
import os


# read a sequence file, ignoring all lines starting with '>' and stripping out non-ACGT characters 
with open("data/sequence/dnaACGT.txt", "r") as file_handle:
    seq = ''.join([filter(lambda x: x in 'acgtACGT', line) for line in file_handle if line[0] != '>'])

# list contents of the matrix directory and check which files look like they should be matrices
matrix_directory = "data/matrix/JASPAR_CORE_2008/"
matrix_names = [filename for filename in os.listdir(matrix_directory) if filename[-4:] == '.pfm']

# MOODS.load_matrix reads JASPAR formated matrix files
matrices = [MOODS.load_matrix(matrix_directory + filename) for filename in matrix_names]

# search with default scoring model
results = MOODS.search(seq, matrices, 0.0001, both_strands=True)

# print the results for each matrix
for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    l = len(matrix[0])
    for (pos, score) in sorted(result, key=lambda (pos,score): pos if pos >= 0 else len(seq) + pos):
        strand = "+"
        hitseq = seq[pos:pos+l]
        if pos < 0:
            pos = len(seq) + pos
            strand = "-"
            hitseq = "".join(reversed(map(lambda s: {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}[s], hitseq.upper())))
        print matrix_name + '|' + str(pos) + '|' + strand + '|' + hitseq + '|'  + str(score)
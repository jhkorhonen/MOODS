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
results = MOODS.search(seq, matrices, 0.001)

# print the results for each matrix
for (matrix_name,results) in zip(matrix_names, results):
    for (pos, score) in results:
        print matrix_name + ':' + str(pos) + ':' + str(score)

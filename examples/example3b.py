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

# search with a custom scoring model
#
# Instead of letting MOODS do log-odds transformation, we do it here 'by hand'.
# What we do here is equivalent to the default MOODS scoring, but other 
# scoring models can be used in the same manner.

pseudocount = 1
p_val = 0.0001

# First, we estimate the bg distribution from the sequence
bg = MOODS.bg_from_sequence(seq, pseudocount)

# Map the matrices to log-odds scoring matrices
matrices = [MOODS.count_log_odds(matrix,bg,pseudocount) for matrix in matrices]

# Thresholds from the p-value
thresholds = [MOODS.threshold_from_p(matrix, bg, p_val) for matrix in matrices]

# Call MOODS search, but tell MOODS not to do log-odds conversion or compute thresholds from a p-value,
# as we have already done all that
results = MOODS.search(seq, matrices, thresholds, convert_log_odds=False, threshold_from_p=False)

# Note that we did not use both_strands=True; as the score depends on the bg distribution, we want to
# score the hits on reverse strand differently. We should rather compute the reverse complement count
# matrices as
#
#  rc_matrices = [reverse_complement(matrix) for matrix in matrices]
#
# and do the log-odds scoring and thresholding for these separately.

# print the results for each matrix
for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
    l = len(matrix[0])
    for (pos, score) in result:
        hitseq = seq[pos:pos+l]
        print matrix_name + '|' + str(pos) + '|'+ hitseq + '|'  + str(score)
import MOODS

import Bio.Seq
import Bio.SeqIO

handle = open("data/sequence/dnaACGT.txt", "r")
records = list(Bio.SeqIO.parse(handle, "fasta"))
handle.close()

seq = records[0]

matrix1 = [     [0,1,0,0,0,0,0,1,1,0],
                [1,0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0,0],
                [0,0,1,1,1,1,1,0,0,1]
            ]
matrix2 = [     [10,0,10,3,5,5],
                [0,5,0,3,5,0,5],
                [0,1,0,3,0,5,0],
                [0,4,0,1,0,0,5]
            ]

results = MOODS.search(seq.seq, [matrix1, matrix2], 0.011)

print("Matrix 1 results: "+ str(len(results[0])))
print("Matrix 2 results: "+ str(len(results[1])))
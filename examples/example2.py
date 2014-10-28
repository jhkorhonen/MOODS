import MOODS

with open("data/sequence/dnaACGT.txt", "r") as file_handle:
    seq = ''.join([filter(lambda x: x in 'acgtACGT', line) for line in file_handle if line[0] != '>'])

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

results = MOODS.search(seq, [matrix1, matrix2], 0.001)

print("Matrix 1 results: "+ str(len(results[0])))
print("Matrix 2 results: "+ str(len(results[1])))
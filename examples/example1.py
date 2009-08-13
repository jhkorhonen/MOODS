import MOODS

matrix = [ [10,0,0],
           [0,10,0],
           [0,0,10],
           [10,10,10]]

results = MOODS.search('actgtggcgtcaacgtaggccaacgtggacccgtacgtaaacgaagaggggtagtc', [matrix], 30, absolute_threshold=30)

for i in results:
    for (position, score) in i:
        print("Position: " + str(position) + " Score: "+ str(score))
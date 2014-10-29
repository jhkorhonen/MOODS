import MOODS

matrix = [ [10,0,0],
           [0,10,0],
           [0,0,10],
           [10,10,10]]

results = MOODS.search('actgtggcgtcaacgtaggccaacgtggacccgtacgtaaacgaagaggggtagtc', [matrix], 30, convert_log_odds=False, threshold_from_p=False)

for i in results:
    for (position, score) in i:
        print("Position: " + str(position) + " Score: "+ str(score))
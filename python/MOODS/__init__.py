import MOODS._cmodule

def loadMatrix(filename):
    array = []
    file = open(filename, "r")
    for line in file:
        array.append(line.split())
    return array

def search(sequence, matrices, thresholds, bg=[0.25, 0.25, 0.25, 0.25], algorithm="lf", q=7, absolute_threshold=False):
    if(type(thresholds) != type(list())):
        lista = []
        for i in matrices:
            lista.append(thresholds)
        thresholds = lista
    return _cmodule.search(str(sequence), matrices, thresholds, bg, algorithm, q, absolute_threshold, True)
import MOODS._cmodule

def loadMatrix(filename):
    '''
    Loads a matrix from a file given.
    Returns an array of arrays of numbers.
    '''
    array = []
    file = open(filename, "r")
    for line in file:
        array.append(line.split())
    return array

def search(sequence, matrices, thresholds, bg=[0.25, 0.25, 0.25, 0.25], algorithm="lf", q=7, absolute_threshold=False):
    '''
    Finds position weight matrix matches in dna sequence. 
    Returns : 
              An array of references to result arrays. There are one result
              array corresponding to each matrix. 
              Each result array is a list of tuples of position and score 
              like: [(pos1, score1), (pos2, score2) ...]
    Parameters: 
             Obligatory:
               sequence  python string object
               matrices
                 One matrix is represented as a typical python multidimensional 
                 array. An array of matrices is expected here.
              thresholds
                 This can be a number, if you want to use the same threshold for
                 all matrices, but you can give multiple thresholds as well. 
                 Multiple thresholds are given as an array of numbers.
             Optional
              bg  
                 Background distribution - an array of four doubles. By default
                 the background is flat.
              algorithm  
                 You can switch search algorithm (Doesn't affect on results)
                    "naive" naive algorithm
                    "pla" permutated lookahead algorithm
                    "supera" super alphabet algorithm. 
                      - Good for long matrices (> 20)
                    "lf" lookahead filtration algorithm. 
                      - Default algorithm in most cases.
                      - Sequence can be searched with multiple matrices 
                        simultaneously. 
                      - You should use this when you have large amount of matrices.
              q  You can optionally give a parameter for algorithm (Doesn't 
                    change results) You can try different values to improve
                    performance.
              absolute_threshold
                    1 if thresholds is given as an absolute value instead of p.
              
    '''
    if(type(thresholds) != type(list())):
        lista = []
        for i in matrices:
            lista.append(thresholds)
        thresholds = lista
    return _cmodule.search(str(sequence), matrices, thresholds, bg, algorithm, q, absolute_threshold, True)
import MOODS._cmodule

def load_matrix(filename):
    '''
    Loads a matrix from a file given.
    Returns matrix as an array of arrays of numbers.
    '''
    array = []
    file = open(filename, "r")
    for line in file:
        array.append(line.split())
    return [map(float, i) for i in array]

def reverse_complement(matrix):
    '''
    Creates a reverse complement of PWM
    '''
    r = [row[:] for row in matrix]
    r.reverse()
    for me in r:
        me.reverse() 
    return r

def transpose(matrix):
    '''
    Creates a transpose of matrix array
    '''
    res = []
    for i in range(len(matrix[0])):
        tmp = []
        for g in matrix:
            tmp.append(g[i])
        res.append(tmp)
    return res

def max_score(matrix):
    '''
    Calculates a maximum score of matrix
    '''
    return sum(map(max, transpose(matrix)))

def bg_from_sequence(seq, ps):
    return _cmodule._bg_from_sequence(str(seq), ps)

def threshold_from_p(matrix, bg, p):
    '''
    Calculates an absolute threshold from a probability value
    '''
    return _cmodule._threshold_from_p(matrix,bg,p)

def count_log_odds(matrix, bg, ps):
    '''
    Calculates a log-odds matrix from a position frequency matrix
    '''
    return _cmodule._count_log_odds(matrix, bg, ps)

def total_matches(matchArray):
    '''
    Calculates a total number of matches.
    '''
    return sum(map((lambda x: len(x)), matchArray))

def flatbg(size = 4):
    '''
    Creates a flat background distribution table
    '''
    return [(1.0/size) for i in range(size)]

def search(sequence, matrices, thresholds, bg=None, algorithm="lf", q=7, absolute_threshold=False, both_strands=False, combine = True):
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
                 the background is calculated from sequence.
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
                    True if thresholds is given as an absolute value instead of p.
              both_strands
                    You can search for both dna strands.
    '''
    if(type(thresholds) != type(list())):
        lista = [thresholds for i in matrices]
        thresholds = lista
    return _cmodule._search(str(sequence), matrices, thresholds, bg, algorithm, q, absolute_threshold, combine, both_strands)

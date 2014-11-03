import MOODS._cmodule
import math

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
    '''
    Estimates the background distribution of nucleotides
    from seq. The pseudocount ps is added to all counts.
    '''
    return _cmodule._bg_from_sequence(str(seq), ps)

def threshold_from_p(matrix, bg, p):
    '''
    Calculates an absolute threshold from a probability value.
    Returns: 
              A threshold value T such that the probability that
              the distribution bg generates a sequence scoring
              at least T on the input matrix is p.
    '''
    return _cmodule._threshold_from_p(matrix,bg,p)

def count_log_odds(matrix, bg, ps, log_base=None):
    '''
    Calculates a log-odds matrix from a position frequency matrix
    Returns: 
              The input PWM matrix transformed to log-odds scores.
              The score for nucleotide N and position i is
              
                 log ( matrix[N][i] + ps * bg[N] ) / C[i]
                   - log ( bg[N] ),
                   
              where C[i] is the sum of terms
              
                 matrix[N][i] + ps * bg[N]
                 
              for all nucleotides N. 
    Parameters: 
             Obligatory:
                matrix
                    Input PWM as a float or integer matrix.
                bg
                    Normalised background distribution given as
                    a list of four floats.
                ps
                    Multiplier for pseudocounts. If matrix is
                    a frequency matrix, you may want to
                    set this to zero.
             Optional:
                log_base  
                    Base for logarithms. Defaults to natural
                    logarithm if None is given.
    '''
    ret = _cmodule._count_log_odds(matrix, bg, ps)
    
    if log_base:
        ret = [[val / math.log(log_base) for val in row] for row in ret]
    
    return ret

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

def search(sequence, matrices, thresholds,
             bg=None,  convert_log_odds=True, threshold_from_p=True, both_strands=False, log_base=None, pseudocount=1,
             algorithm="lf", q=7, combine = True):
    '''
    Finds position weight matrix matches in DNA sequence. 
    Returns: 
              An array of references to result arrays. There is one result
              array for each matrix, in the same order as the input matrices.
              Each result array is a list of tuples of position and score 
              given as: [(pos1, score1), (pos2, score2) ...]
    Parameters: 
             Obligatory:
                sequence
                    DNA sequence as python string object, containing characters
                    acgtACGT.
                matrices
                    An array of matrices, each represented as a list of four
                    lists of equal length. These lists correspond the
                    frequencies or scores of the nucleotides A, C, G and T,
                    respectively.
                thresholds
                    A number or a list of numbers, used as threshold values for
                    matrix scanning.  If a single number is given, it is used
                    for all matrices; otherwise, there should be as many
                    threshold values as there are matrices.
             Optional:
                bg  
                    Background distribution as an array of four doubles,
                    corresponding to the frequencies of A, C, G and T,
                    respectively. By default the background is estimated from
                    the sequence.
                convert_log_odds
                    If True, assumes that the input matrices are frequency or
                    count matrices, and converts them to log-odds scoring
                    matrices using function count_log_odds; otherwise, treat
                    them as scoring matrices. Default True.
                threshold_from_p
                    If True, assumes that thresholds are p-values and computes
                    the corresponding absolute threshold based on the matrix
                    using function threshold_from_p; otherwise the threshold
                    is used as a hard cut-off. Default True.
                log_base
                    Base for logarithms used in log-odds computations. Relevant
                    if using convert_log_odds=True and threshold_from_p=False.
                    Defaults to natural logarithm if None is given.
                pseudocount
                    Pseudocount used in log-odds conversion and added to
                    sequence symbol counts when estimating the background
                    from sequence. Default 1.
                both_strands
                    Scans against reverse complement sequence in addition to
                    the input sequence. Hits on reverse complement are reported
                    at position [position - sequence_length], which is always
                    negative. The actual hit site for any hit is always
                    seq[pos, pos + matrix_length].
                    Default False.
             Tuning parameters:
                (Optional, do not affect the results, but can give minor
                 speed-ups in some cases. You can pretty much ignore these.)
                algorithm  
                     Selects the algorithm to use for scanning
                        "naive" naive algorithm
                        "pla" permutated lookahead algorithm
                        "supera" super alphabet algorithm. 
                          - Good for long matrices (> 20)
                        "lf" lookahead filtration algorithm. 
                          - Default algorithm in most cases.
                          - Sequence can be searched with multiple matrices 
                            simultaneously. 
                q  
                    An integer, used for fine-tuning "supera" and "lf" algorithms.
                    The default value 7 should be ok pretty much always, but can 
                    be tuned to possibly slightly increase performance. 
                combine
                    True or False, determines whether "lf" algorithm combines all
                    matrices to a single scanning pass. Default True.

    '''
    n = len(sequence)
    m = len(matrices)
    
    if not bg:
        bg = bg_from_sequence(sequence, pseudocount)
    elif type(bg) != type(list()) or len(bg) != 4:
        raise RuntimeError('Background does not seem to be a list of length 4')
        return None
    
    if type(thresholds) != type(list()):
        thresholds = [thresholds for i in matrices]
    elif not len(thresholds) == m:
        raise RuntimeError('Number of thresholds does not match the number of matrices')
        return None
    
    if both_strands:
        matrices = matrices + [reverse_complement(matrix) for matrix in matrices]
        thresholds = 2 * thresholds
    
    if convert_log_odds:
        matrices = [count_log_odds(matrix,bg,pseudocount,log_base) for matrix in matrices]
        
    
    if threshold_from_p:
        thresholds = [MOODS.threshold_from_p(matrices[i], bg, thresholds[i]) for i in xrange(len(thresholds))]


    ret = _cmodule._search(str(sequence), matrices, thresholds, bg, algorithm, q, combine)
    
    if both_strands:
        main_hits = ret[:m]
        reverse_hits = [[(pos - n, score) for pos,score in hits] for hits in ret[m:]]
        ret = [main + reverse for main, reverse in zip(main_hits, reverse_hits)]
    
    return ret
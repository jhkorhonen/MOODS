#ifndef MOTIF_H
#define MOTIF_H


#include "moods.h"

namespace MOODS { namespace scan{

class Motif {
    score_matrix mat;
    std::vector<unsigned int> lookahead_order;
    std::vector<double> lookahead_scores;
    
    unsigned int l; // window size
    unsigned int m;
    unsigned int a; // actual alphabet size
    unsigned int av; // 'virtual' alphabet size, first power of 2 >= a
    
    unsigned int wp; // window position
    double T; 
public:
    Motif (const score_matrix& matrix, const vector<double>& bg, unsigned int window_size, double threshold);
    bool window_match(bits_t seq, bits_t shift);
    double window_score(bits_t seq, bits_t shift);
    unsigned int size() { return m; }
    unsigned int alphabet() { return a; }
    unsigned int window_pos() { return l; }
    double threshold() { return T; }
};

}}

#endif
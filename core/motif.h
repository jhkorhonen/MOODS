#ifndef MOTIF_H
#define MOTIF_H


#include "moods.h"

#include <utility>

namespace MOODS { namespace scan{

class Motif {
    score_matrix mat;
    std::vector<unsigned int> lookahead_order;
    std::vector<double> lookahead_scores;
    
    unsigned int l; // window size
    unsigned int m;
    unsigned int a; // actual alphabet size
    
    unsigned int wp; // window position
    double T; 
public:
    Motif (const score_matrix& matrix, const vector<double>& bg, unsigned int window_size, double threshold);
    std:pair<bool, double> window_match(bits_t seq, bits_t shift);
    double Motif::check_hit(const std::vector<unsigned char>& seq, position_t window_match_pos, double score)
    unsigned int size() { return m; }
    unsigned int alphabet() { return a; }
    unsigned int window_pos() { return l; }
    double threshold() { return T; }
};

}}

#endif
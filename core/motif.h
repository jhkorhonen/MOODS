#ifndef MOTIF_H
#define MOTIF_H


#include "moods.h"

#include <utility>

using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{

class Motif {
public:
    virtual std::pair<bool, double> window_match(bits_t seq, bits_t shift) = 0;
    virtual std::pair<bool, double> check_hit(const std::string& s,
                                              const std::vector<unsigned char>& alphabet_map,
                                              const std::size_t window_match_pos,
                                              double score) = 0;
    virtual unsigned int size()=0;
    virtual unsigned int alphabet_size() = 0;
    virtual unsigned int window_pos() = 0;
    virtual double threshold() = 0;
};


// standard 0-order PWM model
class Motif0 : public Motif {
public:
    Motif0 (const score_matrix& matrix, const vector<double>& bg, unsigned int window_size, double threshold);

    std::pair<bool, double> window_match(bits_t seq, bits_t shift);
    std::pair<bool, double> check_hit(const std::string& s, const std::vector<unsigned char>& alphabet_map, const std::size_t window_match_pos, double score);

    unsigned int size() { return m; }
    unsigned int alphabet_size() { return a; }
    unsigned int window_pos() { return wp; }
    double threshold() { return T; }
private:
    score_matrix mat;
    std::vector<unsigned int> lookahead_order;
    std::vector<double> lookahead_scores;
    
    unsigned int l; // window size
    unsigned int m; // length
    unsigned int a; // alphabet size
    
    unsigned int wp; // window position
    double T; 
};

// higher order model with dependencies between adjacent positions
class MotifH : public Motif {
public:
    MotifH (const score_matrix& matrix,
            const vector<double>& bg,
            unsigned int window_size,
            double threshold,
            unsigned int alphabet_size,
            unsigned int dep_len);

    std::pair<bool, double> window_match(bits_t seq, bits_t shift);
    std::pair<bool, double> check_hit(const std::string& s,
                                      const std::vector<unsigned char>& alphabet_map,
                                      const std::size_t window_match_pos,
                                      double score);

    unsigned int size() { return m; }
    unsigned int alphabet_size() { return a; }
    unsigned int window_pos() { return wp; }
    double threshold() { return T; }
private:
    vector<double> expected_scores(const vector<double> &bg);
    vector<double> max_scores_f(size_t start, size_t end);
    vector<double> max_scores_b(size_t start, size_t end);
    vector<double> window_position(const vector<double>& es){
    vector<double> max_prefix_scores();
    vector<double> max_suffix_scores();

    score_matrix mat;
    std::vector<unsigned int> lookahead_order;
    std::vector<double> lookahead_scores;
    
    unsigned int l; // window size
    unsigned int m; // length (of the underlying sequence)
    unsigned int cols; // m - q + 1
    unsigned int rows; 
    unsigned int a; // alphabet size
    unsigned int q; // q-gram length

    vector<double> P; // prefix scores
    vector<double> S; // suffix scores for in-window testing

    
    unsigned int wp; // window position
    double T; 
};

}}

#endif
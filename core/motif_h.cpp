// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include "moods.h"
#include "motif.h"
#include "moods_misc.h"

using std::vector;
using std::size_t;


namespace MOODS { namespace scan{



vector<double> MotifH::expected_scores(const vector<double> &bg){

    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t MASK = (1 << (SHIFT)) - 1;

    vector<double> ret(cols, 0);
    for (int i = 0; i < cols; ++i){
        for (int j = 0; j < rows; ++j){
            double bg_prop = 1;
            for (unsigned int k = 0; k < q; ++k){
                bg_prop *= bg[MASK & (j >> (SHIFT * (q - 1 - k)))];
            }
            ret[i] += bg_prop * mat[j][i];

        }
    }

    return ret;
}

vector<double> MotifH::max_scores_f(size_t start, size_t end){
    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));
    const bits_t Q_MASK = Q_CODE_SIZE - 1;

    vector<double> p_scores(Q_CODE_SIZE, 0);
        
    for (unsigned int i = start; i < end; ++i){
        vector<double> p_scores_n(Q_CODE_SIZE, -std::numeric_limits<double>::infinity());
        
        for (unsigned int j = 0; j < rows; ++j){
            p_scores_n[j & Q_MASK] = std::max(mat[j][i] + p_scores[j >> SHIFT], p_scores_n[j & Q_MASK]);
        }
        p_scores = p_scores_n;
    }

    return p_scores;
}

vector<double> MotifH::max_scores_b(size_t start, size_t end){
    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));
    const bits_t Q_MASK = Q_CODE_SIZE - 1;

    vector<double> p_scores(Q_CODE_SIZE, 0);
        
    for (unsigned int i = 0; i < end-start; ++i){
        vector<double> p_scores_n(Q_CODE_SIZE, -std::numeric_limits<double>::infinity());
        
        for (unsigned int j = 0; j < rows; ++j){
            p_scores_n[j >> SHIFT] = std::max(mat[j][end-i-1] + p_scores[j & Q_MASK], p_scores_n[j >> SHIFT]);
        }
        p_scores = p_scores_n;
    }

    return p_scores;

}

size_t MotifH::window_position(const vector<double>& es){

    if (l >= m){
        return 0;
    }

    // first, we want to compute the max score in each possible window
    // m - l + 1 is the last possible windows position
    vector<double> window_scores(m - l + 1, -std::numeric_limits<double>::infinity());

    for (unsigned int start = 0; start < m - l + 1; ++start){
        vector<double> ss = this->max_scores_f(start, start + l - q - 2);
        
        for (unsigned int j = 0; j < ss.size(); ++j){
            window_scores[start] = std::max(ss[j], window_scores[start]);
        }
    }

    // then we just pick the best window
    double current_exp = 0;
    for (unsigned int i = 0; i < l - q + 1; ++i){
        current_exp = es[i];
    }
    double max_loss = window_scores[0] - current_exp;
    size_t window_pos = 0;
    for (unsigned int i = 1; i < m - l + 1; ++i){
        current_exp -= es[i];
        current_exp += es[i+l-q];

        if (window_scores[i] - current_exp > max_loss){
            window_pos = i+1;
            max_loss = window_scores[i] - current_exp;
        }
    }

    return window_pos;
}

// prefix scores up to the window position
vector<double> MotifH::max_prefix_scores(){
    return this->max_scores_f(0,wp);
}

vector<double> MotifH::max_suffix_scores(){
    return this->max_scores_b(wp+l-q+1,cols);
}




// struct row_comp
// {
//     const vector<double> *ed;
//     bool operator() (int i, int j)
//     {
//         return ( (*ed)[i] > (*ed)[j] );
//     }
// };

// vector<unsigned int> compute_lookahead_order(const vector<double> &ed, unsigned int l, unsigned int window_pos, unsigned int m)
// {
//     if (l >= m)
//     {
//         return vector<unsigned int>();
//     }
//     else
//     {
//         vector<unsigned int> order(m-l, 0);
//         for (int i = 0; i < window_pos; ++i)
//         {
//             order[i] = i;
//         }
//         for (int i = window_pos+l; i < m; ++i)
//         {
//             order[i-l] = i;
//         }
        
//         row_comp comp;
//         comp.ed = &(ed);
        
//         std::sort(order.begin(), order.end(), comp);
        
//         return order;
//     }
// }

// vector<double> compute_lookahead_scores(const score_matrix &mat, const vector<unsigned int> &order, unsigned int l, unsigned int m, unsigned int a)
// {
//     if (l >= m)
//     {
//         return vector<double>();
//     }
//     else
//     {
//         std::vector<double> scores(m-l,0);
        
//         double total = 0;
//         for (int i = m-l-1; i >= 0; --i)
//         {
//             double max = -std::numeric_limits<double>::infinity();
//             for (unsigned int j = 0; j < a; ++j)
//             {
//                 max = std::max(max, mat[j][order[i]]);
//             }
//             total += max;
//             scores[i] = total;
//         }
//         return scores;
//     }
// }

MotifH::MotifH (const score_matrix& matrix, const vector<double>& bg, unsigned int window_size, double threshold, unsigned int alphabet_size)
{

    for (size_t i = 0; i < matrix.size(); ++i){
        for (size_t j = 0; j < matrix[0].size(); ++j){
        }
    }

    mat = matrix;
    T = threshold;

    l = window_size;
    a = alphabet_size;
    cols = mat[0].size();
    rows = mat.size();

    q = MOODS::misc::q_gram_size(rows, a);
    m = mat[0].size() + q - 1;

    vector<double> es = this->expected_scores(bg);
    wp = this->window_position(es);
    P = this->max_prefix_scores();
    S = this->max_suffix_scores();
    
    // lookahead_order = compute_lookahead_order(ed, l, wp, m);
    
    // lookahead_scores = compute_lookahead_scores(mat, lookahead_order, l, m, a);
}

std::pair<bool, double> MotifH::window_match(bits_t seq, bits_t shift)
{
    double score = 0;
    bits_t SHIFT = MOODS::misc::shift(a);
    bits_t MASK = (1 << (SHIFT * q)) - 1;
    
    if (l >= m){
        for (unsigned int i = 0; i < cols; ++i)
        {
            bits_t c = MASK & (seq >> (SHIFT * (l - q - i)));
            score += mat[c][i];
        }
        return std::make_pair(score >= T, score);
    }
    else {
        for (unsigned int i = 0; i < l-q+1; ++i)
        {
            bits_t c = MASK & (seq >> (SHIFT * (l - q - i)));
            score += mat[c][wp+i];
        }

        bits_t prefix = seq >> (SHIFT * (l - q + 1)); // first q - 1 "characters"
        bits_t suffix = seq & ((1 << (SHIFT * (q-1))) - 1); // last q - 1 "characters"
        return std::make_pair(score + P[prefix] + S[suffix] >= T, score);
    }
    
}

std::pair<bool, double> MotifH::check_hit(const std::string& s, const vector<unsigned char>& alphabet_map, const std::size_t window_match_pos, double score)
{
    if (m <= l){
        return  std::make_pair(true, score); // matrix fits fully to the window, so the window score is what we wanted...
    }
    
    // ideally we'd want to do a fancy permuted lookahead like with 0-order models
    // but for now we'll just check the positions in order...

    bits_t SHIFT = MOODS::misc::shift(a);
    bits_t MASK = (1 << (SHIFT * q)) - 1;

    size_t ii = window_match_pos - wp;

    // there's some stuff before the window
    if (wp > 0){
        bits_t CODE = 0;

        for (size_t i = 0; i < q-1; ++i){
            CODE = MASK & ((CODE << SHIFT) ^ alphabet_map[s[ii + i]]);
        }
        for (size_t i = 0; i < wp; ++i){
            CODE = MASK & ((CODE << SHIFT) ^ alphabet_map[s[ii + i + q - 1]]);
            score += mat[CODE][i];
        }
    }

    // stuff after the window
    if (wp < m - l + 1){
        bits_t CODE = 0;

        for (size_t i = wp+l-q+1; i < wp+l; ++i){
            CODE = MASK & (CODE << SHIFT) ^ alphabet_map[s[ii + i]];
        }
        for (size_t i = wp+l-q+1; i < cols; ++i){
            CODE = MASK & (CODE << SHIFT) ^ alphabet_map[s[ii + i + q - 1]];
            score += mat[CODE][i];
        }
    }
    
    return std::make_pair(score >= T, score);
}

} // namespace scan
} // namespace MOODS
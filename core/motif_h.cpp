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

#include <algorithm>

using std::vector;
using std::size_t;

namespace MOODS { namespace scan{



vector<double> MotifH::expected_scores(const vector<double> &bg){

    const bits_t A_MASK = (1 << (SHIFT)) - 1;

    vector<double> ret(cols, 0);
    for (int i = 0; i < cols; ++i){
        for (int j = 0; j < rows; ++j){
            double bg_prop = 1;
            for (unsigned int k = 0; k < q; ++k){
                bg_prop *= bg[A_MASK & (j >> (SHIFT * (q - 1 - k)))];
            }
            ret[i] += bg_prop * mat[j][i];

        }
    }

    return ret;
}

vector<vector<double>> MotifH::max_scores_f(size_t start, size_t end){

    double wsize = end - start;
    vector<vector<double>> max_scores (wsize, vector<double> (Q_CODE_SIZE, 0));

    if (end > start){
        for (unsigned int j = 0; j < rows; ++j){
            max_scores[0][j & Q_MASK] = std::max(mat[j][start], max_scores[0][j & Q_MASK]);
        }

        for (unsigned int i = 1 ; i < wsize; ++i){
            for (unsigned int j = 0; j < rows; ++j){
                max_scores[i][j & Q_MASK] = std::max(mat[j][i+start] + max_scores[i-1][j >> SHIFT], max_scores[i][j & Q_MASK]);
            }
        }
    }

    return max_scores;
}

vector<vector<double>> MotifH::max_scores_b(size_t start, size_t end){

    double wsize = end - start;
    vector<vector<double>> max_scores (wsize, vector<double> (Q_CODE_SIZE, 0));


    vector<double> p_scores(Q_CODE_SIZE, 0);
    if (end > start){

        for (unsigned int j = 0; j < rows; ++j){
            max_scores[wsize-1][j >> SHIFT] = std::max(mat[j][end-1], max_scores[wsize-1][j >> SHIFT]);
        }

        for (unsigned int i = 1; i < wsize; ++i){
            
            for (unsigned int j = 0; j < rows; ++j){
                max_scores[wsize-i-1][j >> SHIFT] = std::max(mat[j][end-i-1] + max_scores[wsize-i][j & Q_MASK], max_scores[wsize-i-1][j >> SHIFT]);
            }
        }
    }

    return max_scores;
}

size_t MotifH::window_position(const vector<double>& es){

    if (l >= m){
        return 0;
    }

    // first, we want to compute the max score in each possible window
    // m - l + 1 is the last possible windows position

    // TODO: fix redundant DP steps
    vector<double> window_scores(m - l + 1, -std::numeric_limits<double>::infinity());

    for (unsigned int start = 0; start < m - l + 1; ++start){
        vector<double> ss = this->max_scores_f(start, start + l - q - 2).back();
        window_scores[start] = *std::max_element(ss.begin(), ss.end());

        // for (unsigned int j = 0; j < ss.size(); ++j){
        //     window_scores[start] = std::max(ss[j], window_scores[start]);
        // }
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
            window_pos = i;
            max_loss = window_scores[i] - current_exp;
        }
    }

    return window_pos;
}

// prefix scores up to the window position
vector<vector<double>> MotifH::max_prefix_scores(){
    return this->max_scores_f(0,wp);
}

vector<vector<double>> MotifH::max_suffix_scores(){
    return this->max_scores_b(wp+l-q+1,cols);
}

MotifH::MotifH (const score_matrix& matrix, const vector<double>& bg, unsigned int window_size, double threshold, unsigned int alphabet_size)
{

    this->mat = matrix;
    this->T = threshold;

    this->l = window_size;
    this->a = alphabet_size;
    this->cols = mat[0].size();
    this->rows = mat.size();

    this->q = MOODS::misc::q_gram_size(rows, a);
    this->m = mat[0].size() + q - 1;

    this->SHIFT = MOODS::misc::shift(a);
    this->MASK = (1 << (SHIFT * q)) - 1;
    this->Q_SHIFT = SHIFT * (q-1);
    this->Q_CODE_SIZE = 1 << Q_SHIFT;
    this->Q_MASK = Q_CODE_SIZE - 1;

    vector<double> es = this->expected_scores(bg);
    this->wp = this->window_position(es);
    this->P = this->max_prefix_scores();
    this->S = this->max_suffix_scores();
    
}

std::pair<bool, double> MotifH::window_match(bits_t seq, bits_t shift)
{
    double score = 0;
    
    if (l >= m){
        for (size_t i = 0; i < cols; ++i)
        {
            bits_t c = MASK & (seq >> (SHIFT * (l - q - i)));
            score += mat[c][i];
        }
        return std::make_pair(score >= T, score);
    }
    else {
        for (size_t i = 0; i < l-q+1; ++i)
        {
            bits_t c = MASK & (seq >> (SHIFT * (l - q - i)));
            score += mat[c][wp+i];
        }

        
        double pot = score;

        if (wp > 0){
            bits_t prefix = seq >> (SHIFT * (l - q + 1)); // first q - 1 "characters"
            pot += P.back()[prefix];
        }
        if (wp < m - l + 1){
            bits_t suffix = seq & ((1 << (SHIFT * (q-1))) - 1); // last q - 1 "characters"
            pot += S.front()[suffix];
        }

        return std::make_pair(pot >= T, score);
    }
    
}

std::pair<bool, double> MotifH::check_hit(const std::string& s, const vector<unsigned char>& alphabet_map, const std::size_t window_match_pos, double score)
{
    if (m <= l){
        return  std::make_pair(true, score); // matrix fits fully to the window, so the window score is what we wanted...
    }
    
    size_t ii = window_match_pos - wp;
    bits_t BACK_CODE = 0;

    // code for the last q-1 positions in the window if we will need that
    if (wp < m - l + 1){
        for (size_t i = 0; i < q-1; ++i){
            BACK_CODE = ( BACK_CODE << SHIFT ) ^ alphabet_map[s[ii + wp + l - q + i + 1]];
        }
    }

    // there's some stuff before the window
    if (wp > 0){
        double forward_threshold = T;

        if (wp < m - l + 1){
            forward_threshold -= S.front()[BACK_CODE];
        }

        bits_t CODE = 0;


        for (size_t i = 0; i < q; ++i){
            CODE = (CODE << SHIFT) ^ alphabet_map[s[ii + wp + i - 1]];
        }

        score += mat[CODE][wp-1];

        // this thing goes backwards
        for (size_t i = 1; i < wp; ++i){

            if (P[wp-i-1][CODE >> SHIFT] + score < forward_threshold){
                return std::make_pair(false, score);
            }

            CODE = (CODE >> SHIFT) ^ (alphabet_map[s[ii + wp - i - 1]] << Q_SHIFT);
            score += mat[CODE][wp-i-1];
        }
    }

    // stuff after the window
    if (wp < m - l + 1){

        // oh look we already precomputed this thing
        bits_t CODE = BACK_CODE;

        // this thing goes forwards
        for (size_t i = wp+l-q+1; i < cols; ++i){

            if (S[i-(wp+l-q+1)][CODE & Q_MASK] + score < T){
                return std::make_pair(false, score);
            }

            CODE = MASK & (CODE << SHIFT) ^ alphabet_map[s[ii + i + q - 1]];
            score += mat[CODE][i];
        }
    }
    return std::make_pair(score >= T, score);
}

} // namespace scan
} // namespace MOODS
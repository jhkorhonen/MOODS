// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include "moods.h"
#include "moods_tools.h"
#include "moods_misc.h"

#include <climits>

using std::vector;
using std::size_t;

namespace MOODS { namespace tools{

// Generates a flat background distribution
vector<double> flat_bg(const unsigned int alphabet_size)
{

    vector<double> bg(alphabet_size, 1.0 / alphabet_size);
    return bg;
}

// Calculates the actual background distribution from a DNA sequence
vector<double> bg_from_sequence_dna(const std::string &seq, const double ps)
{
    unsigned int alphabet_size = 4;
    vector<double> bg(alphabet_size);
    vector<unsigned int> counts(5, 0);
    
    char c;
    // char value;
    
    vector<unsigned char> m(256, 4);
    
    m[(unsigned char)'a'] = 0;
    m[(unsigned char)'A'] = 0;
    
    m[(unsigned char)'c'] = 1;
    m[(unsigned char)'C'] = 1;
    
    m[(unsigned char)'g'] = 2;
    m[(unsigned char)'G'] = 2;
    
    m[(unsigned char)'t'] = 3;
    m[(unsigned char)'T'] = 3;    
    
    for (size_t i = 0; i < seq.size(); ++i){
        c = (unsigned char)seq[i];
        counts[m[c]]++;
    }
    
    for (unsigned int j = 0; j < alphabet_size; ++j){
        bg[j] = (((double)counts[j] + ps)/ ((double)seq.size() + alphabet_size * ps));
    }

    return bg;
}

score_matrix reverse_complement(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    score_matrix ret(a, vector<double>(n));

    for(size_t i = 0; i < a; i++) {
        for (size_t j = 0; j < n; j++) {
            ret[i][j] = mat[a - i - 1][n - j - 1];
        }
    }
    return ret;
}

// Transforms a weight matrix into a PSSM
score_matrix log_odds(const score_matrix &mat, const vector<double> &bg, const double ps)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    score_matrix ret(a, vector<double>(n));

    for (size_t i = 0; i < n; ++i)
    {
        double column_sum = 0;
        for (size_t j = 0; j < a; ++j){
            column_sum += mat[j][i] + ps*bg[j];
        }
        for (size_t j = 0; j < a; ++j) {
            ret[j][i] = log((mat[j][i] + ps*bg[j])/column_sum) - log(bg[j]);
        }
    }
    return ret;
}

// Transforms a weight matrix into a PSSM with non-natural logarithm
score_matrix log_odds(const score_matrix &mat, const vector<double> &bg, const double ps, const double log_base)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    score_matrix ret(a, vector<double>(n));

    for (size_t i = 0; i < n; ++i)
    {
        double column_sum = 0;
        for (size_t j = 0; j < a; ++j){
            column_sum += mat[j][i] + ps*bg[j];
        }
        for (size_t j = 0; j < a; ++j) {
            ret[j][i] = (log((mat[j][i] + ps*bg[j])/column_sum) - log(bg[j])) / log(log_base);
        }
    }
    return ret;
}

// // Calculates a threshold for a scoring matrix from a given p value
double threshold_from_p(const score_matrix &pssm, const vector<double> &bg, const double &p)
{
    const double PVAL_DP_MULTIPLIER = 1000.0;
    
    // Approximate the scoring matrix with integer matrix
    // 'cos we calculate threshold with dynamic programming!
    size_t a = pssm.size();
    size_t n = pssm[0].size();


    vector<vector<int> > mat(a, vector<int>(n));

    int maxT = 0;
    int minV = INT_MAX;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            if (pssm[j][i] > 0.0){
                mat[j][i] = (int) ( PVAL_DP_MULTIPLIER * pssm[j][i] + 0.5 );
            }
            else {
                mat[j][i] = (int) ( PVAL_DP_MULTIPLIER * pssm[j][i] - 0.5 );
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        int max = mat[0][i];
        int min = max;
        for (size_t j = 1; j < a; ++j)
        {
            int v = mat[j][i];
            if (max < v)
                max = v;
            else if (min > v)
                min = v;
        }
        maxT += max;
        if (minV > min)
            minV = min;
    }

    int R = maxT - n * minV;

    vector<double> table0(R + 1, 0.0);
    vector<double> table1(R + 1, 0.0);

    for (size_t j = 0; j < a; ++j)
        table0[mat[j][0] - minV] += bg[j];

    for (size_t i = 1; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            int s = mat[j][i] - minV;
            for (int r = s; r <= R; ++r)
                table1[r] += bg[j] * table0[r - s];
        }
        for (int r = 0; r <= R; ++r)
        {
            table0[r] = table1[r];
            table1[r] = 0.0;
        }
    }

    double sum = 0.0;

    for (int r = R; r >= 0; --r)
    {
        sum += table0[r];
        if (sum > p)
        {
            return (double) ((r + n * minV + 1) / PVAL_DP_MULTIPLIER);
        }
    }

    return (double) ((n * minV) / PVAL_DP_MULTIPLIER);
}


// Calculates maximum possible score for a PSSM
double max_score(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    double ret = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double max = -std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < a; ++j)
        {
            max = std::max(mat[j][i], max);
        }
        ret += max;
    }
    return ret;
}

double min_score(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    double ret = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double min = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < a; ++j)
        {
            min = std::min(mat[j][i], min);
        }
        ret += min;
    }
    return ret;
}

score_matrix log_odds(const vector<vector<double>> &mat, const vector<vector<double>>& low_order_terms,
                      const vector<double> &bg, double ps, size_t a)
{
    size_t rows = mat.size();
    size_t cols = mat[0].size();

    size_t q = MOODS::misc::q_gram_size(rows, a);

    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));

    vector<vector<double>> ret(rows, vector<double>(cols));

    // max-order terms
    for (size_t i = 0; i < cols; ++i){
        for (size_t CODE = 0; CODE < Q_CODE_SIZE; ++CODE){
            double column_sum = 0;
            for (size_t j = 0; j < a; ++j){
                column_sum += mat[(CODE << SHIFT) | j][i] + ps*bg[j];
            }
            for (size_t j = 0; j < a; ++j){
                ret[(CODE << SHIFT) | j][i] = log( (mat[(CODE << SHIFT) | j][i] + ps*bg[j]) / column_sum) - log(bg[j]);
            }
        }        
    }

    // lower-order terms for the first position

    for (size_t r = 0; r < q-1; ++r){
        for (size_t PCODE = 0; PCODE < 1 << (SHIFT * r); ++PCODE ){
            double column_sum = 0;
            for (size_t j = 0; j < a; ++j){
                column_sum += low_order_terms[r][(PCODE << SHIFT) | j] + ps*bg[j];
            }
            for (size_t j = 0; j < a; ++j){
                double lo = log( (low_order_terms[r][(PCODE << SHIFT) | j] + ps*bg[j]) / column_sum) - log(bg[j]);
                size_t PREFIX = ((PCODE << SHIFT) | j) << (SHIFT * (q - r - 1));
                for (size_t SCODE = 0; SCODE < 1 << (SHIFT * (q - r - 1)); ++SCODE){
                    ret[PREFIX | SCODE][0] += lo;
                }
            }
        }
    }

    return ret;
}

score_matrix log_odds_rc(const vector<vector<double>> &mat, const vector<vector<double>>& low_order_terms,
                      const vector<double> &bg, double ps, size_t a){

    size_t q = MOODS::misc::q_gram_size(mat.size(), a);
    size_t rows = mat.size();
    size_t cols = mat[0].size();

    vector<double> bg_rc (bg.size(), 0);

    for (size_t i = 0; i < bg.size(); ++i){
        bg_rc[bg.size() - i - 1] = bg[i];
    }

    score_matrix lo = log_odds(mat, low_order_terms, bg_rc, ps, a);
    score_matrix rc (rows, vector<double> (cols, 0));

    for (size_t i = 0; i < cols; ++i){
        for (size_t j = 0; j < rows; ++j){
            rc[misc::rc_tuple(j, a, q)][cols - i - 1] = lo[j][i];
        }
    }

    return rc;
}


double max_score(const score_matrix &mat, size_t a){

    size_t rows = mat.size();
    size_t cols = mat[0].size();

    unsigned int q = MOODS::misc::q_gram_size(rows, a);

    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));
    const bits_t Q_MASK = Q_CODE_SIZE - 1;

    vector<double> p_scores(Q_CODE_SIZE, 0);
        
    for (unsigned int i = 0; i < cols; ++i){
        vector<double> p_scores_n(Q_CODE_SIZE, -std::numeric_limits<double>::infinity());
        
        for (unsigned int j = 0; j < rows; ++j){
            p_scores_n[j & Q_MASK] = std::max(mat[j][i] + p_scores[j >> SHIFT], p_scores_n[j & Q_MASK]);
        }
        p_scores = p_scores_n;
    }

    double best = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < p_scores.size(); ++i){
        best = std::max(best, p_scores[i]);
    }

    return best;
}

double min_score(const score_matrix &mat, size_t a){

    size_t rows = mat.size();
    size_t cols = mat[0].size();

    unsigned int q = MOODS::misc::q_gram_size(rows, a);

    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));
    const bits_t Q_MASK = Q_CODE_SIZE - 1;

    vector<double> p_scores(Q_CODE_SIZE, 0);
        
    for (unsigned int i = 0; i < cols; ++i){
        vector<double> p_scores_n(Q_CODE_SIZE, std::numeric_limits<double>::infinity());
        
        for (unsigned int j = 0; j < rows; ++j){
            p_scores_n[j & Q_MASK] = std::min(mat[j][i] + p_scores[j >> SHIFT], p_scores_n[j & Q_MASK]);
        }
        p_scores = p_scores_n;
    }

    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < p_scores.size(); ++i){
        best = std::min(best, p_scores[i]);
    }
    
    return best;
}

// temporary threshold-from-p for high-order pwms
double threshold_from_p(const score_matrix &pssm, const vector<double> &bg, const double &p, size_t a)
{
    const double PVAL_DP_MULTIPLIER = 1000.0;
    
    size_t rows = pssm.size();
    size_t cols = pssm[0].size();

    unsigned int q = MOODS::misc::q_gram_size(rows, a);

    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t A_MASK = (1 << SHIFT) - 1;
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));
    const bits_t Q_MASK = Q_CODE_SIZE - 1;



    vector<vector<int> > mat(rows, vector<int>(cols));

    int maxT = 0;
    int minV = INT_MAX;

    for (size_t i = 0; i < cols; ++i)
    {
        for (size_t CODE = 0; CODE < rows; ++CODE)
        {
            if (pssm[CODE][i] > 0.0){
                mat[CODE][i] = (int) ( PVAL_DP_MULTIPLIER * pssm[CODE][i] + 0.5 );
            }
            else {
                mat[CODE][i] = (int) ( PVAL_DP_MULTIPLIER * pssm[CODE][i] - 0.5 );
            }
        }
    }


    for (size_t i = 0; i < cols; ++i)
    {
        int max = mat[0][i];
        int min = max;
        for (size_t CODE = 1; CODE < rows; ++CODE)
        {
            int v = mat[CODE][i];
            if (max < v)
                max = v;
            else if (min > v)
                min = v;
        }
        maxT += max;
        if (minV > min)
            minV = min;
    }

    int R = maxT - cols * minV;

    vector<vector<double>> table0(Q_CODE_SIZE, vector<double>(R + 1, 0));

    for (size_t CODE = 0; CODE < (1 << (SHIFT * q)); ++CODE){
        // TODO: check correctness later for non-DNA alphs
        double prob = 1;

        for (size_t i = 0; i < q; ++i){
            prob *= bg[(CODE >> (q - i - 1)) & A_MASK];
        }

        table0[CODE & Q_MASK][mat[CODE][0] - minV] += prob;
    }

    for (size_t i = 1; i < cols; ++i)
    {
        vector<vector<double>> table1(Q_CODE_SIZE, vector<double>(R + 1, 0.0));
        for (size_t CODE = 0; CODE < (1 << (SHIFT * q)); ++CODE)
        {
            bits_t CODE_PREFIX = (CODE >> SHIFT) & Q_MASK;
            bits_t CODE_SUFFIX = CODE & Q_MASK;
            bits_t CHAR = CODE & A_MASK;
            int s = mat[CODE][i] - minV;
            for (int r = s; r <= R; ++r)
                table1[CODE_SUFFIX][r] += bg[CHAR] * table0[CODE_PREFIX][r - s];
        }
        table0 = table1;
    }


    vector<double> table2(R+1, 0.0);
    for (int r = 0; r < R+1; ++r){
        for (bits_t CODE = 0; CODE < Q_CODE_SIZE; ++CODE){
            table2[r] += table0[CODE][r];
        }
    }


    double sum = 0.0;

    for (int r = R; r >= 0; --r)
    {
        sum += table2[r];
        if (sum > p)
        {
            return (double) ((r + cols * minV + 1) / PVAL_DP_MULTIPLIER);
        }
    }

    return (double) ((cols * minV) / PVAL_DP_MULTIPLIER);
}


} // namespace tools
} // namespace MOODS
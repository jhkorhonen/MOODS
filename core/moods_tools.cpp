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
#include "match_types.h"

#include <limits>

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
    
    unsigned int count_all = counts[0] + counts[1] + counts[2] + counts[3];
    
    for (unsigned int j = 0; j < alphabet_size; ++j){
        bg[j] = (((double)counts[j] + ps)/ ((double)count_all + alphabet_size * ps));
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
    
    score_matrix ret = log_odds(mat, bg, ps);
    
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j) {
            ret[j][i] = ret[j][i] / log(log_base);
        }
    }
    
    return ret;

    
    //
    // score_matrix ret(a, vector<double>(n));
    //
    // for (size_t i = 0; i < n; ++i)
    // {
    //     double column_sum = 0;
    //     for (size_t j = 0; j < a; ++j){
    //         column_sum += mat[j][i] + ps*bg[j];
    //     }
    //     for (size_t j = 0; j < a; ++j) {
    //         ret[j][i] = (log((mat[j][i] + ps*bg[j])/column_sum) - log(bg[j])) / log(log_base);
    //     }
    // }
    // return ret;
}

// // Calculates a threshold for a scoring matrix from a given p value
double threshold_from_p_with_precision(const score_matrix &pssm, const vector<double> &bg, const double &p, const double precision)
{
    
    // Approximate the scoring matrix with integer matrix for DP
    long a = pssm.size();
    long n = pssm[0].size();


    vector<vector<long> > mat(a, vector<long>(n));

    long maxT = 0;
    long minV = std::numeric_limits<long>::max();

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            if (pssm[j][i] > 0.0){
                mat[j][i] = (long) ( precision * pssm[j][i] + 0.5 );
            }
            else {
                mat[j][i] = (long) ( precision * pssm[j][i] - 0.5 );
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        long max = mat[0][i];
        long min = max;
        for (size_t j = 1; j < a; ++j)
        {
            long v = mat[j][i];
            if (max < v)
                max = v;
            else if (min > v)
                min = v;
        }
        maxT += max;
        if (minV > min)
            minV = min;
    }

    long R = maxT - n * minV;

    vector<double> table0(R + 1, 0.0);
    vector<double> table1(R + 1, 0.0);

    for (size_t j = 0; j < a; ++j)
        table0[mat[j][0] - minV] += bg[j];

    for (size_t i = 1; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            long s = mat[j][i] - minV;
            for (long r = s; r <= R; ++r)
                table1[r] += bg[j] * table0[r - s];
        }
        for (long r = 0; r <= R; ++r)
        {
            table0[r] = table1[r];
            table1[r] = 0.0;
        }
    }


        
    double sum = table0[R];

    if (sum > p){
        return max_score(pssm) - min_delta(pssm)/2;
    }
        
    for (long r = R-1; r >= 0; --r)
    {
        sum += table0[r];
        if (sum > p)
        {
            return (double) ((r + n * minV + 1) / precision);
        }
    }

    return min_score(pssm)-1.0;
}

double threshold_from_p(const score_matrix &pssm, const vector<double> &bg, const double &p)
{
    return threshold_from_p_with_precision(pssm, bg, p, DEFAULT_DP_PRECISION);
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

double min_delta(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    double min_delta = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i)
    {
        double max = -std::numeric_limits<double>::infinity();
        double second_max = -std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < a; ++j)
        {
            if (mat[j][i] > max){
                second_max = max;
                max = mat[j][i];
            }
            else if (mat[j][i] < max && mat[j][i] > second_max) {
                second_max = mat[j][i];
            }
        }
        min_delta = std::min(min_delta, max - second_max);
    }
    return min_delta;
}

score_matrix log_odds(const vector<vector<double>> &mat, const vector<vector<double>>& low_order_terms,
                      const vector<double> &bg, double ps, const size_t a)
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

score_matrix log_odds(const vector<vector<double>> &mat, const vector<vector<double>>& low_order_terms,
                      const vector<double> &bg, double ps, const size_t a, const double log_base)
{
    size_t rows = mat.size();
    size_t cols = mat[0].size();
    
    score_matrix ret = log_odds(mat, low_order_terms, bg, ps, a);
    
    for (size_t i = 0; i < cols; ++i)
    {
        for (size_t j = 0; j < rows; ++j) {
            ret[j][i] = ret[j][i] / log(log_base);
        }
    }
    
    return ret;
}


score_matrix reverse_complement(const vector<vector<double>> &mat, const size_t a){

    size_t q = MOODS::misc::q_gram_size(mat.size(), a);
    size_t rows = mat.size();
    size_t cols = mat[0].size();
    
    score_matrix ret(rows, vector<double>(cols));
    
    for (size_t i = 0; i < cols; ++i){
        for (size_t j = 0; j < rows; ++j){
            ret[misc::rc_tuple(j, a, q)][cols - i - 1] = mat[j][i];
        }
    }

    return ret;
}


double max_score(const score_matrix &mat, const size_t a){

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

double min_score(const score_matrix &mat, const size_t a){

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
double threshold_from_p_with_precision(const score_matrix &pssm, const vector<double> &bg, const double &p, const double precision, const size_t a)
{
    
    long rows = pssm.size();
    long cols = pssm[0].size();

    unsigned int q = MOODS::misc::q_gram_size(rows, a);

    const bits_t SHIFT = MOODS::misc::shift(a);
    const bits_t A_MASK = (1 << SHIFT) - 1;
    const bits_t Q_CODE_SIZE  =  (1 << (SHIFT * (q-1)));
    const bits_t Q_MASK = Q_CODE_SIZE - 1;

    vector<vector<long> > mat(rows, vector<long>(cols));

    long maxT = 0;
    long minV = std::numeric_limits<long>::max();

    for (size_t i = 0; i < cols; ++i)
    {
        for (size_t CODE = 0; CODE < rows; ++CODE)
        {
            if (pssm[CODE][i] > 0.0){
                mat[CODE][i] = (long) ( precision * pssm[CODE][i] + 0.5 );
            }
            else {
                mat[CODE][i] = (long) ( precision * pssm[CODE][i] - 0.5 );
            }
        }
    }


    for (size_t i = 0; i < cols; ++i)
    {
        long max = mat[0][i];
        long min = max;
        for (size_t CODE = 1; CODE < rows; ++CODE)
        {
            long v = mat[CODE][i];
            if (max < v)
                max = v;
            else if (min > v)
                min = v;
        }
        maxT += max;
        if (minV > min)
            minV = min;
    }

    long R = maxT - cols * minV;

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
            long s = mat[CODE][i] - minV;
            for (long r = s; r <= R; ++r)
                table1[CODE_SUFFIX][r] += bg[CHAR] * table0[CODE_PREFIX][r - s];
        }
        table0 = table1;
    }


    vector<double> table2(R+1, 0.0);
    for (long r = 0; r < R+1; ++r){
        for (bits_t CODE = 0; CODE < Q_CODE_SIZE; ++CODE){
            table2[r] += table0[CODE][r];
        }
    }


    double sum = table2[R];
    
    if (sum > p){
        return max_score(pssm, a) - min_delta(pssm)/2;;
    }

    for (long r = R-1; r >= 0; --r)
    {
        sum += table2[r];
        if (sum > p)
        {
            return (double) ((r + cols * minV + 1) / precision);
        }
    }

    return min_score(pssm, a) - 1.0;
}

double threshold_from_p(const score_matrix &pssm, const vector<double> &bg, const double &p, const size_t a)
{
    return threshold_from_p_with_precision(pssm, bg, p, DEFAULT_DP_PRECISION, a);
}


vector<MOODS::variant> snp_variants(const std::string &seq){
    
    vector<MOODS::variant> ret;    
    vector<std::string> snp_alt(256);
    
    snp_alt[(unsigned char)'w'] = "AT";
    snp_alt[(unsigned char)'W'] = "AT";
    
    snp_alt[(unsigned char)'s'] = "CG";
    snp_alt[(unsigned char)'S'] = "CG";
    
    snp_alt[(unsigned char)'m'] = "AC";
    snp_alt[(unsigned char)'M'] = "AC";

    snp_alt[(unsigned char)'k'] = "GT";    
    snp_alt[(unsigned char)'K'] = "GT";

    snp_alt[(unsigned char)'r'] = "AG";    
    snp_alt[(unsigned char)'R'] = "AG";
    
    snp_alt[(unsigned char)'y'] = "CT";
    snp_alt[(unsigned char)'Y'] = "CT";
    
    snp_alt[(unsigned char)'b'] = "CGT";
    snp_alt[(unsigned char)'B'] = "CGT";
    
    snp_alt[(unsigned char)'d'] = "AGT";
    snp_alt[(unsigned char)'D'] = "AGT";
    
    snp_alt[(unsigned char)'h'] = "ACT";
    snp_alt[(unsigned char)'H'] = "ACT";
    
    snp_alt[(unsigned char)'v'] = "ACG";
    snp_alt[(unsigned char)'V'] = "ACG";
    
    for (size_t i = 0; i < seq.size(); ++i){
        for (size_t j = 0; j < snp_alt[seq[i]].size(); ++j){
            ret.push_back(variant{i,i+1,snp_alt[seq[i]].substr(j,1)});
        }
    }
    
    return ret;
}



} // namespace tools
} // namespace MOODS
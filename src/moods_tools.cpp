// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include <vector>
#include <cmath>

using std::vector;

namespace MOODS { namespace tools{

// Generates a flat background distribution
vector<double> flat_bg(const unsigned int alphabet_size)
{

    vector<double> bg(alphabet_size, 1.0 / alphabet_size);
    return bg;
}

// Calculates the actual background distribution from the sequence
vector<double> bg_from_sequence(const vector<unsigned char> &seq, const int alphabet_size, const double ps)
{
    vector<double> bg(alphabet_size);
    vector<unsigned int> counts(alphabet_size, 0);

    for (unsigned int i = 0; i < seq.size(); ++i){
        ++counts[seq[i]];
    }
    for (int j = 0; j < size; ++j){
        bg[j] = (((double)counts[j] + ps)/ ((double)seq.size() + alphabet_size * ps));
    }
    
    return bg;
}

score_matrix reverse_complement(score_matrix &mat) {
    
    unsigned int a = mat.size();
    unsigned int n = mat[0].size();
    
    score_matrix ret(a, vector<score_t>(n));
    
    for(unsigned int i = 0; i < a; i++) {
        for (unsigned int j = 0; j < n; j++) {
            ret[i][j] = mat[a - i - 1][n - j - 1];
        }
    }
    return ret;
}

// Transforms a weight matrix into a PSSM
score_matrix log_odds(const scoreMatrix &mat, const vector<double> &bg, const double ps)
{
    int a = mat.size();
    int n = mat[0].size();

    score_matrix ret(a, vector<score_t>(n));

    for (int i = 0; i < n; ++i)
    {
        double column_sum = 0;
        for (int j = 0; j < a; ++j){
            column_sum += mat[j][i] + ps*bg[j];
        }
        for (int j = 0; j < a; ++j) {
            ret[j][i] = log((mat[j][i] + ps*bg[j])/column_sum) - log(bg[j]);
        }
    }
    return ret;
}

// Transforms a weight matrix into a PSSM with non-natural logarithm
scoreMatrix log_odds(const score_matrix &mat, const vector<double> &bg, const double ps, const double log_base)
{
    int a = mat.size();
    int n = mat[0].size();

    score_matrix ret(a, vector<score_t>(n));

    for (int i = 0; i < n; ++i)
    {
        double column_sum = 0;
        for (int j = 0; j < a; ++j){
            column_sum += mat[j][i] + ps*bg[j];
        }
        for (int j = 0; j < a; ++j) {
            ret[j][i] = (log((mat[j][i] + ps*bg[j])/column_sum) - log(bg[j])) / log(log_base);
        }
    }
    return ret;
}

// Calculates a threshold for a scoring matrix from a given p value
score_t tresholdFromP(const scoreMatrix &pssm, const doubleArray &bg, const double &p)
{
    // Approximate the scoring matrix with integer matrix
    // 'cos we calculate threshold with dynamic programming!
    int numA = pssm.size();
    int n = pssm[0].size();

    intArray col(n, 0);
    intMatrix mat(numA, col);

    int maxT = 0;
    int minV = INT_MAX;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < numA; ++j)
        {
            if (pssm[j][i] > 0.0){
                mat[j][i] = (int) ( PVAL_DP_MULTIPLIER * pssm[j][i] + 0.5 );
            }
            else {
                mat[j][i] = (int) ( PVAL_DP_MULTIPLIER * pssm[j][i] - 0.5 );
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        int max = mat[0][i];
        int min = max;
        for (int j = 1; j < numA; ++j)
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

    doubleArray table0(R + 1, 0.0);
    doubleArray table1(R + 1, 0.0);

    for (int j = 0; j < numA; ++j)
        table0[mat[j][0] - minV] += bg[j];

    for (int i = 1; i < n; ++i)
    {
        for (int j = 0; j < numA; ++j)
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
            return (score_t) ((r + n * minV + 1) / PVAL_DP_MULTIPLIER);
        }
    }

	return (score_t) ((n * minV) / PVAL_DP_MULTIPLIER);
}

// Calculates maximum possible score for a PSSM
score_t max_score(const score_matrix &mat)
{
    int a = mat.size();
    int n = mat[0].size();
    
    score_t ret = 0;
    for (int i = 0; i < n; ++i)
    {
        score_t max = SCORE_MIN;
        for (int j = 0; j < a; ++j)
        {
            max = std::max(mat[j][i], max);
        }
        ret += max;
    }
    return ret;
}

score_t min_score(const score_matrix &mat)
{
    int a = mat.size();
    int n = mat[0].size();
    
    score_t ret = 0;
    for (int i = 0; i < n; ++i)
    {
        score_t max = SCORE_MAX;
        for (int j = 0; j < a; ++j)
        {
            max = std::min(mat[j][i], max);
        }
        ret += max;
    }
    return ret;
}



} // namespace tools
} // namespace MOODS
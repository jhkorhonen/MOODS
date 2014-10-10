// pssm_algorithms - collection of algorithms for finding PSSM matches from sequences
// Copyright (C) 2007-2009  Pasi Rastas, Janne Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.



#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>
#include <list>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <cassert>
#include <limits.h>

using std::cin;
using std::cout;
using std::cerr;
using std::vector;
using std::list;
using std::queue;
using std::deque;


#include "pssm_algorithms.hpp"

// *************************
// Data structures and types
// *************************

// Struct for comparing matrix rows with sort()
struct compareRows
{
    const doubleArray *goodness;
    bool operator() (int i, int j)
    {
        return ( (*goodness)[i] > (*goodness)[j] );
    }
};

// Output list element for mm AC automaton
struct OutputListElementMulti
{
    score_t score;
    int matrix;
    bool full;
};

// ****************
// Helper functions
// ****************

// Generates a flat background distribution
doubleArray flatBG(const int size)
{

    doubleArray bg(size, 0);
    for (int i = 0; i < size; ++i)
    {
        bg[i] = 1.0 / size;
    }
    return bg;
}

doubleArray bgFromSequence(SeqIterator &seq_it, const int size, const double pseudocount)
{
    doubleArray bg(size);
    intArray weights(size, 0);

    for (int i = 0; i < size; ++i)
        weights[i] = 0;

    unsigned long seq_size = 0;
    while(seq_it.hasData()) {
    	++weights[*seq_it];
    	++seq_it;
    	seq_size ++;
    }


    for (int i = 0; i < size; ++i)
    {
        bg[i] = (((double)weights[i] + pseudocount)/ ((double)seq_size + size * pseudocount));
    }

    return bg;
}

// Calculates the actual background distribution from the sequence
doubleArray bgFromSequence(const charArray &seq, const int size, const double pseudocount)
{
    doubleArray bg(size);
    intArray weights(size, 0);

    for (int i = 0; i < size; ++i)
        weights[i] = 0;
    for (unsigned int i = 0; i < seq.size(); ++i)
        ++weights[seq[i]];

    for (int i = 0; i < size; ++i)
    {
        bg[i] = (((double)weights[i] + pseudocount)/ ((double)seq.size() + size * pseudocount));
    }

    return bg;
}

charArray readString(std::basic_istream<char> &in)
{
    charArray ret;
    while ( !in.eof() )
    {
        int value = 0;
        in >> value;
        ret.push_back((char) value);
    }
    return ret;
}

// Reads a dna sequence from input
charArray readDNA(std::basic_istream<char> &in)
{
    charArray ret;
    while ( !in.eof() )
    {
        char c = 0;
        in >> c;
        int value = 0;
        switch (c)
        {
            case 'a':
            case 'A': value = 0; break;
            case 'c':
            case 'C': value = 1; break;
            case 'g':
            case 'G': value = 2; break;
            case 't':
            case 'T': value = 3; break;
            case 'n':
            case 'N': value = (int) (4 * rand() / (1.0 + RAND_MAX)); break;
            case '>':   // Skip FASTA header lines
                while (c != '\n' && !in.eof())
                {
                    c = in.get();
                }
                value = 5;
                break;
            default:
                value = 5;
                break;
        }
        if (value <= 3)
            ret.push_back((char) value);
    }
    return ret;
}


// Reads a character sequence from input; accepts letters A-Z
charArray readAZ(std::basic_istream<char> &in)
{
    charArray ret;
    while ( !in.eof() )
    {
        char c = 0;
        in >> c;
        int value = c - 'A';
        if (value >= 0 && value < 26)
            ret.push_back((char) value);
    }
    return ret;
}

// Reads a matrix from input
scoreMatrix readMatrix(std::basic_istream<char> &in)
{
    scoreMatrix ret;
    unsigned int row = 0;
    char ch = 0;
    bool read = false;
    while ( in.get(ch) )
    {
        if (ch == '\n')
        {
            if (read)
            {
                ++row;
                read = false;
            }
        }
        else
        {
            in.putback(ch);
            score_t value = 0;
            in >> value;
            if (row == ret.size())
                ret.push_back(scoreArray());
            ret[row].push_back(value);
            read = true;
        }
    }
    return ret;
}

scoreMatrix reverseComplement(scoreMatrix &mat) {
	if(mat.size() != 4)
		return scoreMatrix();

	scoreMatrix ret;

	for(int i = 3; i >= 0; i--) {
		scoreArray row;
		for(int k = mat[i].size()-1; k>=0; k--) {
			row.push_back(mat[i][k]);
		}
		ret.push_back(row);
	}
	return ret;
}

// Transforms a weight matrix into a PSSM
scoreMatrix counts2LogOdds(const scoreMatrix &mat, const doubleArray &bg, const double ps)
{
    int numA = mat.size();
    int n = mat[0].size();

    vector<score_t> col(n, 0);
    scoreMatrix ret(numA, col);

    for (int i = 0; i < n; ++i)
    {

        int count = 0;
        for (int j = 0; j < numA; ++j)
            count += mat[j][i];

        for (int j = 0; j < numA; ++j)
        {
            double f = (mat[j][i] + bg[j])  / (count + ps);
            ret[j][i] = log(f) - log(bg[j]);
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
score_t maxScore(const scoreMatrix &mat)
{
    int numA = mat.size();
    int n = mat[0].size();
    score_t ret = 0;
    for (int i = 0; i < n; ++i)
    {
        score_t max = SCORE_MIN;
        for (int j = 0; j < numA; ++j)
        {
            score_t ma = mat[j][i];
            if (ma > max)
                max = ma;
        }
        ret += max;
    }
    return ret;
}

doubleArray expectedDifferences(const scoreMatrix &mat, const doubleArray &bg)
{
    int numA = mat.size();
    int m = mat[0].size();
    doubleArray ret(m);

    for (int i = 0; i < m; ++i)
    {
        score_t max = SCORE_MIN;
        for (int j = 0; j < numA; ++j)
        {
            if (max < mat[j][i])
                max = mat[j][i];
        }

        ret[i] = max;

        for (int j = 0; j < numA; ++j)
        {
            ret[i] -= bg[j] * mat[j][i];
        }
    }

    return ret;
}


void printScoreMatrix(const scoreMatrix &mat)
{
    int numA = mat.size();
    int n = mat[0].size();

    for (int j = 0; j < numA; ++j)
    {
        for (int i = 0; i < n; ++i)
            cout << std::setw(8) << std::fixed << std::setprecision(3) << mat[j][i];
        cout << "\n";
    }
}

void printSequence(const charArray &seq, const int n)
{
    for (int i = 0; i < n; ++i)
        cout << (int) seq[i] << " ";
    cout << "\n";
}

void printMatchArray(const matchArray &match)
{
    for (unsigned int i = 0; i < match.size(); ++i)
    {
        cout << match[i].position << "\t" << match[i].score << "\n";
    }

}

// *******************
// Matching algorithms
// *******************


matchArray naiveAlgorithm(SeqIterator &s, const scoreMatrix &p, const score_t tol)
{
    const int m = p[0].size();
    matchArray ret;
    matchData hit;
    for(;s.hasData(m); ++s) {
    	score_t score = 0;
		for(int j= 0; j < m; ++j)
		{
			score += p[ s[j] ][j];
		}
		if (score >= tol)
		{
			hit.position = s.position();
			hit.score = score;
			ret.push_back(hit);
		}
    }
    return ret;
}

// Permutated lookahead algorithm
matchArray permutatedLookAhead(SeqIterator &seq_it, const scoreMatrix &p, const doubleArray &bg, const score_t tol)
{
    const int numA = p.size();
    const int m = p[0].size();

    intArray order(m, 0);

    doubleArray row_goodness = expectedDifferences(p, bg);


    for(unsigned int i = 0; i < order.size(); ++i)
    {
        order[i] = i;
    }

    compareRows comp;
    comp.goodness = &row_goodness;

    sort(order.begin(), order.end(), comp);

    scoreArray T(m, 0);
    for (int j = m - 1; j > 0; --j)
    {
	    score_t max = SCORE_MIN;
	    for (int i = 0; i < numA; ++i)
	    {
	        if (max < p[i][order[j]])
                max = p[i][order[j]];
	    }
	    T[j - 1] = T[j] + max;
    }

    intArray::iterator x;
    intArray::iterator start = order.begin();

    matchArray ret;
    matchData hit;


    score_t score = 0;
    for (; seq_it.hasData(m); ++seq_it)
    {
        score = 0;
        x = start;
        for (int j = 0; j < m; ++j)
        {
            score = score + p[ seq_it[*x] ][(*x)];
            if (score + T[j] < tol)
                break;
            ++x;
        }
        if (score >= tol)
        {
            hit.position = seq_it.position();
            hit.score = score;
            ret.push_back(hit);
        }
    }
    return ret;
}

// Naive superalphabet algorithm - good for very long matrices
// for DNA
matchArray naiveSuperalphabetAlgorithmDNA(const int q, SeqIterator &seq_it, const scoreMatrix &p, const score_t tol)
{

    const int m = p[0].size();


    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

   // const position_t n = s.size();

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;

    const int mod = m % q; // counts length of end part...
    int numQ2 = m / q;
    if (mod != 0)
	   ++numQ2;

    const int numQ = numQ2;

    bitArray sA(q, 0);
    bits_t code = 0;

    scoreMatrix scores;

    for (int i = 0; i < numQ; ++i)
	   scores.push_back(scoreArray(size, 0));

    {
	int k = 0;
	for (int i = 0; i < m; i+=q)
	{
            score_t score = 0;
            for (int j = 0; j < q; ++j)
                if ( (i + j) < m )
                    score = score + p[0][i + j];
            scores[k++][0] = score;
        }
    }

    while (true)
    {
        ++code;
        if (code >= size)
            break;

        for (int i = 0; i < numQ; ++i)
            scores[i][code] = scores[i][code - 1];

        int pos = q - 1;
        while (pos >= 0)
        {
            if (sA[pos] < numA - 1)
            {
                bits_t j = 0;
                bits_t sa = sA[pos];
                for (int i = 0; i < m - mod; i+=q)
                    scores[j++][code] += (p[sa + 1][pos + i] - p[sa][pos + i]);
                if (mod != 0 && pos < mod)
                {
                    int i = m - mod;
                    scores[j][code] += (p[sa + 1][pos + i] - p[sa][pos + i]);
                }
                ++sA[pos];
                break;
            }
            else
            {
                bits_t j = 0;
                for (int i = 0; i < m - mod; i+=q)
                    scores[j++][code] += (p[0][pos + i] - p[numA - 1][pos + i]);
                if (mod != 0 && pos < mod)
                {
                    int i = m - mod;
                    scores[j][code] += (p[0][pos + i] - p[numA - 1][pos + i]);
                }
                sA[pos] = 0;
                --pos;
            }
        }
    }

    bitArray codes(numQ, 0);

    if(!seq_it.hasData(q + m)) {
    	return matchArray();
    }
    int k = 0;
    for (int i = 0; i < m; i+=q)
    {
        bits_t c = 0;
        bits_t digit = (q - 1) * BITSHIFT;

        for (int ii = 0; ii < q - 1; ++ii)
        {
            digit -= BITSHIFT;
            c += (seq_it[i + ii] << digit);
        }
        codes[k++] = c;
    }

    matchArray ret;
    matchData hit;

	position_t pos;
	score_t score;
	bits_t c;
    for ( ;seq_it.hasData(numQ*q); ++seq_it)
    {
        pos = q - 1;

        score = 0;
        for (int j = 0; (j < numQ); ++j)
        {
            c = codes[j];
            c = ((c << BITSHIFT) + (seq_it[pos])) & BITAND;
            score += scores[j][c];
            codes[j] = c;
            pos += q;
        }
        if (score >= tol)
        {
            hit.position = seq_it.position();
            hit.score = score;
            ret.push_back(hit);
        }
    }
    // A hack for handling last few positions of the sequence
	for (; seq_it.hasData(m); ++seq_it)
	{
		score_t score = 0;
		for (int j = 0; j < m; ++j)
		{
			score += p[ seq_it[j]][j];
		}
		if (score >= tol)
		{
			hit.position = seq_it.position();
			hit.score = score;
			ret.push_back(hit);
		}
	}

    return ret;
}

// LFA for DNA
matchArray lookaheadFiltrationDNA(const int q, SeqIterator &seq_it, const scoreMatrix &p, const doubleArray &bg, const score_t tol)
{
    const int m = p[0].size();

	if (q > m)
    {
        return lookaheadFiltrationDNA(m, seq_it, p, bg, tol);
    }
    const int BITSHIFT = 2;
    const bits_t numA = 4; // 2**BITSIFT
    //const position_t n = s.size();

    // Find the window

    double current_goodness = 0;
    doubleArray goodness = expectedDifferences(p, bg);

    for (int i = 0; i < q; ++i)
    {
        current_goodness += goodness[i];
    }

    double max_goodness = current_goodness;
    int window_pos = 0;

    for (int i = 0; i < m - q; ++i)
    {
        current_goodness -= goodness[i];
        current_goodness += goodness[i+q];
        if (current_goodness > max_goodness)
        {
            max_goodness = current_goodness;
            window_pos = i+1;
        }
    }
    // Arrange matrix indeces not in window by entropy, for use in scanning
    intArray order(m-q, 0);
    for (int i = 0; i < window_pos; ++i)
    {
        order[i] = i;
    }
    for (int i = window_pos+q; i < m; ++i)
    {
        order[i-q] = i;
    }
    compareRows comp;
    comp.goodness = &goodness;
    sort(order.begin(), order.end(), comp);


	// Lookahead array for indeces outside the window
    scoreArray L(m-q+1, 0);
    for (int j = m-q; j > 0; --j)
    {
	    score_t max = SCORE_MIN;
	    for (int i = 0; i < (int) numA; ++i)
	    {
	        if (max < p[i][order[j-1]])
                max = p[i][order[j-1]];
	    }
	    L[j - 1] = L[j] + max;
    }

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;

    bitArray sA(q, 0);
    bits_t code = 0;


	// Scores in the window for all possible strings
    scoreArray scores(size, 0);

    {
    score_t score = 0;
    for (int j = 0; j < q; ++j)
        score = score + p[0][window_pos + j];
    scores[0] = score;
    }

    while (true)
    {
        ++code;
        if (code >= size)
            break;

        scores[code] = scores[code - 1];

        int pos = q - 1;
        while (pos >= 0)
        {
            if (sA[pos] < numA - 1)
            {
                bits_t sa = sA[pos];
                scores[code] += (p[sa + 1][pos+window_pos] - p[sa][pos+window_pos]);
                ++sA[pos];
                break;
            }
            else
            {
                scores[code] += (p[0][pos+window_pos] - p[numA - 1][pos+window_pos]);
                sA[pos] = 0;
                --pos;
            }
        }
    }

	// Actual scanning

	matchArray ret;
	matchData hit;

    intArray::iterator y = order.begin();
    position_t offset = q + window_pos - 1;

    score_t limit = tol - L[0];
    score_t score = 0;

    if(!seq_it.hasData(offset)) {
    	assert(false);
    	return ret;
    }

    code = 0;
    for (position_t ii = window_pos; ii < offset; ++ii)
        code = (code << BITSHIFT) + seq_it[ii];

    for (; seq_it.hasData(m); ++seq_it)
    {
        code = ((code << BITSHIFT) + seq_it[offset]) & BITAND;
        if (scores[code] >= limit)
        {
            score = scores[code];
            y = order.begin();
            for (int j = 0; j < m-q  ;++j)
            {
                if (score + L[j] < tol)
                    break;
                score += p[seq_it[(*y)]][*y];
                ++y;
            }
            if (score >= tol){
				hit.position = seq_it.position();
				hit.score = score;
				ret.push_back(hit);
			}
        }
    }

	return ret;

}

// Multimatrix LFA for DNA
vector<matchArray> multipleMatrixLookaheadFiltrationDNA(const int q, SeqIterator &seq_it, const vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol)
{
    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

    intArray m(matrices.size(), 0);
    for (int i = 0; i < (int) matrices.size(); ++i)
    {
        m[i] = matrices[i][0].size();
    }

    // Calculate entropies for all matrices
    vector<doubleArray> goodnesses;
    goodnesses.reserve(matrices.size());

    for (int i = 0; i < (int)matrices.size(); ++i)
    {
        goodnesses.push_back(expectedDifferences(matrices[i], bg));
    }

    intArray window_positions;
    window_positions.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            window_positions.push_back(0);
        }
        else
        {
            double current_goodness = 0;
            for (int i = 0; i < q; ++i)
            {
                current_goodness += goodnesses[k][i];
            }

            double max_goodness = current_goodness;
            int window_pos = 0;

            for (int i = 0; i < m[k] - q; ++i)
            {
                current_goodness -= goodnesses[k][i];
                current_goodness += goodnesses[k][i+q];
                if (current_goodness > max_goodness)
                {
                    max_goodness = current_goodness;
                    window_pos = i+1;
                }
            }
            window_positions.push_back(window_pos);
        }
    }

    // Calculate lookahead scores for all matrices
    scoreMatrix T;
    T.reserve(matrices.size());

    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        scoreArray C(m[k],0);
        for (int j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T.push_back(C);
    }

    // Pre-window scores
    scoreArray P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        score_t B = 0;
        for (int j = 0; j < window_positions[k]; ++j)
        {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            B += max;
        }
        P.push_back(B);
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    intMatrix orders;
    orders.reserve(matrices.size());
    scoreMatrix L;
    L.reserve(matrices.size());
    for (unsigned short k = 0; k < (int) matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            intArray temp1;
            orders.push_back(temp1);
            scoreArray temp2;
            L.push_back(temp2);
        }
        else
        {
            intArray order(m[k]-q, 0);
            for (int i = 0; i < window_positions[k]; ++i)
            {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i)
            {
                order[i-q] = i;
            }

            compareRows comp;
            comp.goodness = &(goodnesses[k]);

            sort(order.begin(), order.end(), comp);

            orders.push_back(order);

            scoreArray K(m[k]-q, 0);
            for (int j = m[k]-q-1; j > 0; --j)
            {
                score_t max = INT_MIN;
                for (unsigned int i = 0; i < numA; ++i)
                {
                    if (max < matrices[k][i][order[j]])
                        max = matrices[k][i][order[j]];
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;


    vector<vector< OutputListElementMulti> > output(size);

	{
		bitArray sA(q,0);
		while (true)
		{
			bits_t code = 0;
			for (int j = 0; j < q; ++j)
			{
				code = (code << BITSHIFT) | sA[j];
			}

			for (unsigned int k = 0; k < matrices.size(); ++k )
			{
                if (m[k] <= q)
                {
                    score_t score = 0;
                    for (int i = 0; i < m[k]; ++i)
                    {
                        score += matrices[k][sA[i]][i];
                    }
                    if (score >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = true;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
                else
                {
                    score_t score = 0;
                    for (int i = 0; i < q; ++i)
                    {
                        score += matrices[k][sA[i]][i + window_positions[k]];
                    }
                    if (score + P[k] + T[k][q + window_positions[k]-1] >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = false;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
            }

			int pos = 0;
			while (pos < q)
			{
				if (sA[pos] < numA - 1)
				{
					++sA[pos];
					break;
				}
				else
				{
					sA[pos] = 0;
					++pos;
				}
			}

			if (pos == q)
				break;
		}
	}

	// Scanning
	vector<matchArray> ret;
    for (unsigned int i = 0; i < matrices.size(); ++i)
    {
        matchArray temp;
        ret.push_back(temp);
    }

	matchData hit;

    score_t score;
    position_t k;
    position_t limit;
    position_t ii;
    score_t tolerance;
    intArray::iterator z;

    bits_t code = 0;
    if(!seq_it.hasData(q - 1)) {
    	assert(false);
    	return ret;
    }

    for (position_t ii = 0; ii < q - 1; ++ii)
        code = (code << BITSHIFT) + seq_it[ii];


	position_t offset = 0;
	position_t max_offset = 0;
	for(position_t kk = 0; kk < (position_t) window_positions.size();  kk++) {
		if(max_offset < window_positions[kk])
			max_offset = window_positions[kk];
	}

	while(seq_it.hasData(q - 1))
	{
		code = ((code << BITSHIFT) + seq_it[offset + q - 1]) & BITAND;
		if (!output[code].empty())
		{
			for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
			{
				if (y->full) // A Hit for a matrix of length <= q
				{
					hit.position = seq_it.position() + offset;
					hit.score = y->score;
					ret[y->matrix].push_back(hit);
					continue;
				}
				if (seq_it.hasData(offset + m[y->matrix] - window_positions[y->matrix])) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
				{
					score = y->score;
					k = y->matrix;
					limit = m[k] - q;
					ii = offset - window_positions[k];
					if(ii < 0)
						continue;

					tolerance = tol[k];
					z = orders[k].begin();
					for (int j = 0; j < limit  ;++j)
					{
						score += matrices[k][seq_it[ii+(*z)]][*z];
						if (score + L[k][j] < tolerance)
							break;
						++z;
					}
					if (score >= tolerance){
						hit.position = seq_it.position() + offset - window_positions[k];
						hit.score = score;
						ret[k].push_back(hit);
					}
				}
			}
		}
		if(offset < max_offset)
			++offset;
		else
			++seq_it;
	}

    seq_it += offset;
	for (; seq_it.hasData(); ++seq_it) // possible hits for matrices shorter than q near the end of sequence
	{
		code = (code << BITSHIFT) & BITAND; // dummy character to the end of code

		if (!output[code].empty())
        {

			for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
        	{
            	if (y->full && seq_it.hasData(m[y->matrix])) // only sufficiently short hits are considered
            	{
                	hit.position = seq_it.position();
                	hit.score = y->score;
                	ret[y->matrix].push_back(hit);
            	}
        	}
		}
	}
	return ret;
}

vector<matchArray> searchDNA(SeqIterator &seq, const vector<scoreMatrix> &matrices, const scoreArray &tol, const doubleArray &bg, const int q) {
	if(matrices.size() <= 0)
		return vector<matchArray>();
	int alphabet_size = matrices[0].size();
	if(alphabet_size != 4) {
		vector<matchArray> matches;
		for(unsigned int i = 0; i < matrices.size(); i++) {
			seq.reset();
			matches.push_back(naiveAlgorithm(seq, matrices[i], tol[i]));
		}
		return matches;
	}
	seq.reset();
	if(matrices.size() == 1) {
		vector<matchArray> matches;
		if(matrices[0].size() > 20) {
			matches.push_back(naiveSuperalphabetAlgorithmDNA(q,seq, matrices[0], tol[0]));
		}
		else {
			matches.push_back(lookaheadFiltrationDNA(q, seq, matrices[0], bg, tol[0]));
		}
		return matches;
	}

	return multipleMatrixLookaheadFiltrationDNA(q, seq, matrices, bg, tol);
}

// *******************
// Matching algorithms char array versions
// *******************

// The naive algorithm
matchArray naiveAlgorithm(const charArray &s, const scoreMatrix &p, const score_t tol)
{
    const int m = p[0].size();
    const position_t n = s.size();

    matchArray ret;
    matchData hit;

    for (position_t i = 0; i <= n - m; ++i)
    {
        score_t score = 0;
        for (int j = 0; j < m; ++j)
        {
            score += p[ s[i + j] ][j];
        }
        if (score >= tol)
        {
            hit.position = i;
            hit.score = score;
            ret.push_back(hit);
        }
    }

    return ret;
}


// Permutated lookahead algorithm
matchArray permutatedLookAhead(const charArray &s, const scoreMatrix &p, const doubleArray &bg, const score_t tol)
{
    const int numA = p.size();
    const int m = p[0].size();
    const position_t n = s.size();

    intArray order(m, 0);

    doubleArray row_goodness = expectedDifferences(p, bg);


    for(unsigned int i = 0; i < order.size(); ++i)
    {
        order[i] = i;
    }

    compareRows comp;
    comp.goodness = &row_goodness;

    sort(order.begin(), order.end(), comp);

    scoreArray T(m, 0);
    for (int j = m - 1; j > 0; --j)
    {
	    score_t max = SCORE_MIN;
	    for (int i = 0; i < numA; ++i)
	    {
	        if (max < p[i][order[j]])
                max = p[i][order[j]];
	    }
	    T[j - 1] = T[j] + max;
    }

    intArray::iterator x;
    intArray::iterator start = order.begin();

    matchArray ret;
    matchData hit;

    charArray::const_iterator y = s.begin();


    score_t score = 0;
    for (position_t i = 0; i <= n - m; ++i)
    {
        score = 0;
        x = start;
        for (int j = 0; j < m; ++j)
        {
            score = score + p[ y[*x] ][(*x)];
            if (score + T[j] < tol)
                break;
            ++x;
        }
        if (score >= tol)
        {
            hit.position = i;
            hit.score = score;
            ret.push_back(hit);
        }
        ++y;
    }

    return ret;
}



// Naive superalphabet algorithm - good for very long matrices
// For arbitary alphabet
matchArray naiveSuperalphabetAlgorithm(const int q, const charArray &s, const scoreMatrix &p, const score_t tol)
{

    const int m = p[0].size();
    if (q > m)
	{
		return naiveSuperalphabetAlgorithm(m, s, p, tol);
    }



    const unsigned int numA = p.size(); // 2**BITSIFT

    const position_t n = s.size();

    int size2 = 1;
    for (int i = 0; i < q; ++i)
		size2 *= numA;
    const bits_t size = size2;

    const int mod = m % q; // counts length of end part...
    int numQ2 = m / q;
    if (mod != 0)
		++numQ2;

    const int numQ = numQ2;

    bitArray sA(q, 0);
    bits_t code = 0;

    scoreMatrix scores;

    for (int i = 0; i < numQ; ++i)
		scores.push_back(scoreArray(size, 0));

	int k = 0;
	for (int i = 0; i < m; i+=q)
	{
            score_t score = 0;
            for (int j = 0; j < q; ++j)
                if ( (i + j) < m )
                    score = score + p[0][i + j];
            scores[k++][0] = score;
    }

    while (true)
	{
		++code;
		if (code >= size)
	    	break;

		for (int i = 0; i < numQ; ++i)
	    	scores[i][code] = scores[i][code - 1];

		int pos = q - 1;
		while (pos >= 0)
		{
	    	if (sA[pos] < numA - 1)
			{
				int j = 0;
				int sa = sA[pos];
				for (int i = 0; i < m - mod; i+=q)
		    		scores[j++][code] += (p[sa + 1][pos + i] - p[sa][pos + i]);
				if (mod != 0 && pos < mod)
				{
		    		int i = m - mod;
		    		scores[j][code] += (p[sa + 1][pos + i] - p[sa][pos + i]);
				}
				++sA[pos];
				break;
	    	}
	    	else
			{
				int j = 0;
				for (int i = 0; i < m - mod; i+=q)
		    		scores[j++][code] += (p[0][pos + i] - p[numA - 1][pos + i]);
				if (mod != 0 && pos < mod)
				{
		    		int i = m - mod;
		    		scores[j][code] += (p[0][pos + i] - p[numA - 1][pos + i]);
				}
				sA[pos] = 0;
				--pos;
	    	}
		}
    }

    bitArray codes(numQ, 0);

    k = 0;
    for (int i = 0; i < m; i+=q)
	{
		bits_t c = 0;
		for (int ii = 0; ii < q - 1; ++ii)
	    	c = (s[i + ii]  +  c * numA);
		codes[k++] = c;
    }

    charArray old(numQ, 0);

	matchArray ret;
	matchData hit;

	position_t pos;
	score_t score;
	bits_t c;

    for (position_t i = 0; i <= n - (numQ*q); ++i)
	{
		pos = i + q - 1;

		score = 0;
		for (int j = 0; j < numQ; ++j)
		{
	    	c = codes[j];
	    	c = (s[pos]  + c * numA - size * old[j]);
	    	old[j] = s[pos - q + 1];
	    	score += scores[j][c];
	    	codes[j] = c;
	    	pos += q;
		}
		if (score >= tol)
		{
			hit.position = i;
            hit.score = score;
            ret.push_back(hit);
		}
    }

    // A hack for handling last few positions of the sequence
    for (position_t i = n - (numQ*q) + 1; i <= n - m; ++i)
    {
        score_t score = 0;
        for (int j = 0; j < m; ++j)
        {
            score += p[ s[i + j] ][j];
        }
        if (score >= tol)
        {
            hit.position = i;
            hit.score = score;
            ret.push_back(hit);
        }
    }

	return ret;
}


// Naive superalphabet algorithm - good for very long matrices
// for DNA
matchArray naiveSuperalphabetAlgorithmDNA(const int q, const charArray &s, const scoreMatrix &p, const score_t tol)
{

    const int m = p[0].size();
    if (q > m)
    {
        return naiveSuperalphabetAlgorithmDNA(m, s, p, tol);
    }

    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

    const position_t n = s.size();

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;

    const int mod = m % q; // counts length of end part...
    int numQ2 = m / q;
    if (mod != 0)
	   ++numQ2;

    const int numQ = numQ2;

    bitArray sA(q, 0);
    bits_t code = 0;

    scoreMatrix scores;

    for (int i = 0; i < numQ; ++i)
	   scores.push_back(scoreArray(size, 0));

    {
	int k = 0;
	for (int i = 0; i < m; i+=q)
	{
            score_t score = 0;
            for (int j = 0; j < q; ++j)
                if ( (i + j) < m )
                    score = score + p[0][i + j];
            scores[k++][0] = score;
        }
    }

    while (true)
    {
        ++code;
        if (code >= size)
            break;

        for (int i = 0; i < numQ; ++i)
            scores[i][code] = scores[i][code - 1];

        int pos = q - 1;
        while (pos >= 0)
        {
            if (sA[pos] < numA - 1)
            {
                bits_t j = 0;
                bits_t sa = sA[pos];
                for (int i = 0; i < m - mod; i+=q)
                    scores[j++][code] += (p[sa + 1][pos + i] - p[sa][pos + i]);
                if (mod != 0 && pos < mod)
                {
                    int i = m - mod;
                    scores[j][code] += (p[sa + 1][pos + i] - p[sa][pos + i]);
                }
                ++sA[pos];
                break;
            }
            else
            {
                bits_t j = 0;
                for (int i = 0; i < m - mod; i+=q)
                    scores[j++][code] += (p[0][pos + i] - p[numA - 1][pos + i]);
                if (mod != 0 && pos < mod)
                {
                    int i = m - mod;
                    scores[j][code] += (p[0][pos + i] - p[numA - 1][pos + i]);
                }
                sA[pos] = 0;
                --pos;
            }
        }
    }

    bitArray codes(numQ, 0);

    int k = 0;
    for (int i = 0; i < m; i+=q)
    {
        bits_t c = 0;
        bits_t digit = (q - 1) * BITSHIFT;

        for (int ii = 0; ii < q - 1; ++ii)
        {
            digit -= BITSHIFT;
            c += (s[i + ii] << digit);
        }
        codes[k++] = c;
    }

    matchArray ret;
    matchData hit;

	position_t pos;
	score_t score;
	bits_t c;

    for (position_t i = 0; i <= n - (numQ*q); ++i)
    {
        pos = i + q - 1;

        score = 0;
        for (int j = 0; j < numQ; ++j)
        {
            c = codes[j];
            c = ((c << BITSHIFT) + (s[pos])) & BITAND;
            score += scores[j][c];
            codes[j] = c;
            pos += q;
        }
        if (score >= tol)
        {
            hit.position = i;
            hit.score = score;
            ret.push_back(hit);
        }
    }

    // A hack for handling last few positions of the sequence
    for (position_t i = n - (numQ*q) + 1; i <= n - m; ++i)
    {
        score_t score = 0;
        for (int j = 0; j < m; ++j)
        {
            score += p[ s[i + j] ][j];
        }
        if (score >= tol)
        {
            hit.position = i;
            hit.score = score;
            ret.push_back(hit);
        }
    }

    return ret;
}


// LFA for arbitary alphabet
matchArray lookaheadFiltration(const int q, const charArray &s, const scoreMatrix &p, const doubleArray &bg, const score_t tol)
{
    const int m = p[0].size();
    if (q > m)
    {
        return lookaheadFiltration(m, s, p, bg, tol);
    }

    const unsigned int numA = p.size();
    const position_t n = s.size();

    // Find the window

    double current_goodness = 0;

    doubleArray goodness = expectedDifferences(p, bg);
    for (int i = 0; i < q; ++i)
    {
        current_goodness += goodness[i];
    }

    double max_goodness = current_goodness;
    int window_pos = 0;

    for (int i = 0; i < m - q; ++i)
    {
        current_goodness -= goodness[i];
        current_goodness += goodness[i+q];
        if (current_goodness > max_goodness)
        {
            max_goodness = current_goodness;
            window_pos = i+1;
        }
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    intArray order(m-q, 0);
    for (int i = 0; i < window_pos; ++i)
    {
        order[i] = i;
    }
    for (int i = window_pos+q; i < m; ++i)
    {
        order[i-q] = i;
    }

    compareRows comp;
    comp.goodness = &goodness;

    sort(order.begin(), order.end(), comp);

    scoreArray L(m-q+1, 0);
    for (int j = m-q; j > 0; --j)
    {
	    score_t max = SCORE_MIN;
	    for (int i = 0; i < (int) numA; ++i)
	    {
	        if (max < p[i][order[j-1]])
                max = p[i][order[j-1]];
	    }
	    L[j - 1] = L[j] + max;
    }

    const bits_t size = (int) pow(numA, q); // numA^q
	const bits_t divisor = (int) pow(numA, q-1);

    bitArray sA(q, 0);
    bits_t code = 0;

    scoreArray scores(size, 0);


    {
    score_t score = 0;
    for (int j = 0; j < q; ++j)
        score = score + p[0][window_pos + j];
    scores[0] = score;
    }

    while (true)
    {
        ++code;
        if (code >= size)
            break;

        scores[code] = scores[code - 1];

        int pos = q - 1;
        while (pos >= 0)
        {
            if (sA[pos] < numA - 1)
            {
                bits_t sa = sA[pos];
                scores[code] += (p[sa + 1][pos+window_pos] - p[sa][pos+window_pos]);
                ++sA[pos];
                break;
            }
            else
            {
                scores[code] += (p[0][pos+window_pos] - p[numA - 1][pos+window_pos]);
                sA[pos] = 0;
                --pos;
            }
        }
    }

    bitArray shifts(size, 0);

	for (bits_t i = 0; i < size; ++i){
		shifts[i] = (i * numA) - ((i / divisor) * size);
	}

    // Actual scanning

	matchArray ret;
	matchData hit;

    intArray::iterator y = order.begin();
    int offset = q + window_pos - 1;

    score_t limit = tol - L[0];
    score_t score = 0;

    code = 0;
    for (position_t ii = window_pos; ii < offset; ++ii)
        code = (code * numA) + s[ii];

    for (position_t i = 0; i <= n - m + 1; ++i)
    {
        code = shifts[code] + s[i + offset];

        if (scores[code] >= limit)
        {
            score = scores[code];
            y = order.begin();
            for (position_t j = 0; j < m-q  ;++j)
            {
                if (score + L[j] < tol)
                    break;
                score += p[s[i+(*y)]][*y];
                ++y;
            }
            if (score >= tol){
				hit.position = i;
				hit.score = score;
				ret.push_back(hit);
			}
        }
    }

    return ret;
}


// LFA for DNA
matchArray lookaheadFiltrationDNA(const int q, const charArray &s, const scoreMatrix &p, const doubleArray &bg, const score_t tol)
{
    const int m = p[0].size();
    if (q > m)
    {
        return lookaheadFiltrationDNA(m, s, p, bg, tol);
    }

    const int BITSHIFT = 2;
    const bits_t numA = 4; // 2**BITSIFT
    const position_t n = s.size();

    // Find the window

    double current_goodness = 0;
    doubleArray goodness = expectedDifferences(p, bg);


    for (int i = 0; i < q; ++i)
    {
        current_goodness += goodness[i];
    }

    double max_goodness = current_goodness;
    int window_pos = 0;

    for (int i = 0; i < m - q; ++i)
    {
        current_goodness -= goodness[i];
        current_goodness += goodness[i+q];
        if (current_goodness > max_goodness)
        {
            max_goodness = current_goodness;
            window_pos = i+1;
        }
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    intArray order(m-q, 0);
    for (int i = 0; i < window_pos; ++i)
    {
        order[i] = i;
    }
    for (int i = window_pos+q; i < m; ++i)
    {
        order[i-q] = i;
    }

    compareRows comp;
    comp.goodness = &goodness;

    sort(order.begin(), order.end(), comp);


	// Lookahead array for indeces outside the window
    scoreArray L(m-q+1, 0);
    for (int j = m-q; j > 0; --j)
    {
	    score_t max = SCORE_MIN;
	    for (int i = 0; i < (int) numA; ++i)
	    {
	        if (max < p[i][order[j-1]])
                max = p[i][order[j-1]];
	    }
	    L[j - 1] = L[j] + max;
    }

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;

    bitArray sA(q, 0);
    bits_t code = 0;


	// Scores in the window for all possible strings
    scoreArray scores(size, 0);

    {
    score_t score = 0;
    for (int j = 0; j < q; ++j)
        score = score + p[0][window_pos + j];
    scores[0] = score;
    }

    while (true)
    {
        ++code;
        if (code >= size)
            break;

        scores[code] = scores[code - 1];

        int pos = q - 1;
        while (pos >= 0)
        {
            if (sA[pos] < numA - 1)
            {
                bits_t sa = sA[pos];
                scores[code] += (p[sa + 1][pos+window_pos] - p[sa][pos+window_pos]);
                ++sA[pos];
                break;
            }
            else
            {
                scores[code] += (p[0][pos+window_pos] - p[numA - 1][pos+window_pos]);
                sA[pos] = 0;
                --pos;
            }
        }
    }

	// Actual scanning

	matchArray ret;
	matchData hit;

    intArray::iterator y = order.begin();
    position_t offset = q + window_pos - 1;

    score_t limit = tol - L[0];
    score_t score = 0;

    code = 0;
    for (position_t ii = window_pos; ii < offset; ++ii)
        code = (code << BITSHIFT) + s[ii];


    for (position_t i = 0; i <= n - m + 1; ++i)
    {
        code = ((code << BITSHIFT) + s[i + offset]) & BITAND;

        if (scores[code] >= limit)
        {
            score = scores[code];
            y = order.begin();
            for (int j = 0; j < m-q  ;++j)
            {
                if (score + L[j] < tol)
                    break;
                score += p[s[i+(*y)]][*y];
                ++y;
            }
            if (score >= tol){
				hit.position = i;
				hit.score = score;
				ret.push_back(hit);
			}
        }
    }

	return ret;

}


// Multmatrix LFA for arbitary alphabet
vector<matchArray> multipleMatrixLookaheadFiltration(const int q, const charArray &s, const vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol)
{

    const unsigned int numA = matrices[0].size();

    const position_t n = s.size();

    intArray m(matrices.size(), 0);
    for (int i = 0; i < (int) matrices.size(); ++i)
    {
        m[i] = matrices[i][0].size();
    }

    // Calculate entropies for all matrices
    vector<doubleArray> goodnesses;
    goodnesses.reserve(matrices.size());

    for (int i = 0; i < (int)matrices.size(); ++i)
    {
        goodnesses.push_back(expectedDifferences(matrices[i], bg));
    }

    intArray window_positions;
    window_positions.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            window_positions.push_back(0);
        }
        else
        {
            double current_goodness = 0;
            for (int i = 0; i < q; ++i)
            {
                current_goodness += goodnesses[k][i];
            }

            double max_goodness = current_goodness;
            int window_pos = 0;

            for (int i = 0; i < m[k] - q; ++i)
            {
                current_goodness -= goodnesses[k][i];
                current_goodness += goodnesses[k][i+q];
                if (current_goodness > max_goodness)
                {
                    max_goodness = current_goodness;
                    window_pos = i+1;
                }
            }
            window_positions.push_back(window_pos);
        }
    }

    // Calculate lookahead scores for all matrices
    scoreMatrix T;
    T.reserve(matrices.size());

    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        scoreArray C(m[k],0);
        for (int j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T.push_back(C);
    }

    // Pre-window scores
    scoreArray P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        score_t B = 0;
        for (int j = 0; j < window_positions[k]; ++j)
        {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            B += max;
        }
        P.push_back(B);
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    intMatrix orders;
    orders.reserve(matrices.size());
    scoreMatrix L;
    L.reserve(matrices.size());

    for (int k = 0; k < (int) matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            intArray temp1;
            orders.push_back(temp1);
            scoreArray temp2;
            L.push_back(temp2);
        }
        else
        {
            intArray order(m[k]-q, 0);
            for (int i = 0; i < window_positions[k]; ++i)
            {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i)
            {
                order[i-q] = i;
            }

            compareRows comp;
            comp.goodness = &(goodnesses[k]);

            sort(order.begin(), order.end(), comp);

            orders.push_back(order);

            scoreArray K(m[k]-q, 0);
            for (int j = m[k]-q-1; j > 0; --j)
            {
                score_t max = SCORE_MIN;
                for (unsigned int i = 0; i < numA; ++i)
                {
                    if (max < matrices[k][i][order[j]])
                        max = matrices[k][i][order[j]];
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }


    // Note that this is not strictly speaking a bit parallel algorithm
    // Assuming bits_t is a unsigned integer type, this doesn't matter
    const bits_t size = (int) pow(numA, q);

    vector<vector< OutputListElementMulti> > output(size);

	{
		bitArray sA(q,0);
		while (true)
		{
			bits_t code = 0;
			for (int j = 0; j < q; ++j)
			{
				code = (code * numA) + sA[j];
			}

			for (unsigned int k = 0; k < matrices.size(); ++k )
			{
                if (m[k] <= q)
                {
                    score_t score = 0;
                    for (int i = 0; i < m[k]; ++i)
                    {
                        score += matrices[k][sA[i]][i];
                    }
                    if (score >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = true;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
                else
                {
                    score_t score = 0;
                    for (int i = 0; i < q; ++i)
                    {
                        score += matrices[k][sA[i]][i + window_positions[k]];
                    }
                    if (score + P[k] + T[k][q + window_positions[k]-1] >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = false;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
            }

			int pos = 0;
			while (pos < q)
			{
				if (sA[pos] < numA - 1)
				{
					++sA[pos];
					break;
				}
				else
				{
					sA[pos] = 0;
					++pos;
				}
			}

			if (pos == q)
				break;
		}
	}

    score_t score;
    int k;
    position_t limit;
    position_t ii;
    score_t tolerance;
    intArray::iterator z;

	vector<matchArray> ret;
    for (unsigned int i = 0; i < matrices.size(); ++i)
    {
        matchArray temp;
        ret.push_back(temp);
    }
	matchData hit;

    bits_t code = 0;
    for (position_t ii = 0; ii < q - 1; ++ii)
        code = (code * numA) + s[ii];

    bits_t last = 0;

    for (position_t i = 0; i < n - q + 1; ++i)
    {
        code = ((code * numA) + s[i + q - 1]) - (size * last);
        last = s[i];

        if (!output[code].empty())
        {
            for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
            {
                if (y->full)
                {
                    hit.position = i;
                    hit.score = y->score;
                    ret[y->matrix].push_back(hit);
                    continue;
                }
                if (i - window_positions[y->matrix] >= 0 && i + m[y->matrix] - window_positions[y->matrix] <= n) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                {
                    score = y->score;
                    k = y->matrix;
                    limit = m[k] - q;
                    ii = i - window_positions[k];
                    tolerance = tol[k];
                    z = orders[k].begin();
                    for (int j = 0; j < limit  ;++j)
                    {
                        score += matrices[k][s[ii+(*z)]][*z];
                        if (score + L[k][j] < tolerance)
                            break;
                        ++z;
                    }
                    if (score >= tolerance){
						hit.position = i - window_positions[k];
						hit.score = score;
						ret[k].push_back(hit);
					}
                }
            }
        }
    }

	for (position_t i = n - q + 1; i < n; ++i) // possible hits for matrices shorter than q near the end of sequence
	{
        code = (code * numA) - (size * last);
        last = s[i];

		if (!output[code].empty())
        {

			for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
        	{
            	if (y->full && m[y->matrix] < n - i + 1) // only sufficiently short hits are considered
            	{
                	hit.position = i;
                	hit.score = y->score;
                	ret[y->matrix].push_back(hit);
            	}
        	}
		}
	}

	return ret;
}



// Multimatrix LFA for DNA
vector<matchArray> multipleMatrixLookaheadFiltrationDNA(const int q, const charArray &s, const vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol)
{

    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

    const position_t n = s.size();

    intArray m(matrices.size(), 0);
    for (int i = 0; i < (int) matrices.size(); ++i)
    {
        m[i] = matrices[i][0].size();
    }

    // Calculate entropies for all matrices
    vector<doubleArray> goodnesses;
    goodnesses.reserve(matrices.size());

    for (int i = 0; i < (int)matrices.size(); ++i)
    {
        goodnesses.push_back(expectedDifferences(matrices[i], bg));
    }

    intArray window_positions;
    window_positions.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            window_positions.push_back(0);
        }
        else
        {
            double current_goodness = 0;
            for (int i = 0; i < q; ++i)
            {
                current_goodness += goodnesses[k][i];
            }

            double max_goodness = current_goodness;
            int window_pos = 0;

            for (int i = 0; i < m[k] - q; ++i)
            {
                current_goodness -= goodnesses[k][i];
                current_goodness += goodnesses[k][i+q];
                if (current_goodness > max_goodness)
                {
                    max_goodness = current_goodness;
                    window_pos = i+1;
                }
            }
            window_positions.push_back(window_pos);
        }
    }

    // Calculate lookahead scores for all matrices
    scoreMatrix T;
    T.reserve(matrices.size());

    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        scoreArray C(m[k],0);
        for (int j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T.push_back(C);
    }

    // Pre-window scores
    scoreArray P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        score_t B = 0;
        for (int j = 0; j < window_positions[k]; ++j)
        {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            B += max;
        }
        P.push_back(B);
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    intMatrix orders;
    orders.reserve(matrices.size());
    scoreMatrix L;
    L.reserve(matrices.size());

    for (unsigned short k = 0; k < (int) matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            intArray temp1;
            orders.push_back(temp1);
            scoreArray temp2;
            L.push_back(temp2);
        }
        else
        {
            intArray order(m[k]-q, 0);
            for (int i = 0; i < window_positions[k]; ++i)
            {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i)
            {
                order[i-q] = i;
            }

            compareRows comp;
            comp.goodness = &(goodnesses[k]);

            sort(order.begin(), order.end(), comp);

            orders.push_back(order);

            scoreArray K(m[k]-q, 0);
            for (int j = m[k]-q-1; j > 0; --j)
            {
                score_t max = INT_MIN;
                for (unsigned int i = 0; i < numA; ++i)
                {
                    if (max < matrices[k][i][order[j]])
                        max = matrices[k][i][order[j]];
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;


    vector<vector< OutputListElementMulti> > output(size);

	{
		bitArray sA(q,0);
		while (true)
		{
			bits_t code = 0;
			for (int j = 0; j < q; ++j)
			{
				code = (code << BITSHIFT) | sA[j];
			}

			for (unsigned int k = 0; k < matrices.size(); ++k )
			{
                if (m[k] <= q)
                {
                    score_t score = 0;
                    for (int i = 0; i < m[k]; ++i)
                    {
                        score += matrices[k][sA[i]][i];
                    }
                    if (score >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = true;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
                else
                {
                    score_t score = 0;
                    for (int i = 0; i < q; ++i)
                    {
                        score += matrices[k][sA[i]][i + window_positions[k]];
                    }
                    if (score + P[k] + T[k][q + window_positions[k]-1] >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = false;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
            }

			int pos = 0;
			while (pos < q)
			{
				if (sA[pos] < numA - 1)
				{
					++sA[pos];
					break;
				}
				else
				{
					sA[pos] = 0;
					++pos;
				}
			}

			if (pos == q)
				break;
		}
	}

	// Scanning

	vector<matchArray> ret;
    for (unsigned int i = 0; i < matrices.size(); ++i)
    {
        matchArray temp;
        ret.push_back(temp);
    }

	matchData hit;

    score_t score;
    position_t k;
    position_t limit;
    position_t ii;
    score_t tolerance;
    intArray::iterator z;

    bits_t code = 0;
    for (position_t ii = 0; ii < q - 1; ++ii)
        code = (code << BITSHIFT) + s[ii];


    for (position_t i = 0; i < n - q + 1; ++i)
    {
    	code = ((code << BITSHIFT) + s[i + q - 1]) & BITAND;

        if (!output[code].empty())
        {
            for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
            {
                if (y->full) // A Hit for a matrix of length <= q
                {
                    hit.position = i;
                    hit.score = y->score;
                    ret[y->matrix].push_back(hit);
                    continue;
                }
                if (i - window_positions[y->matrix] >= 0 && i + m[y->matrix] - window_positions[y->matrix] <= n) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                {
                    score = y->score;
                    k = y->matrix;
                    limit = m[k] - q;
                    ii = i - window_positions[k];
                    tolerance = tol[k];
                    z = orders[k].begin();
                    for (int j = 0; j < limit  ;++j)
                    {
                        score += matrices[k][s[ii+(*z)]][*z];
                        if (score + L[k][j] < tolerance)
                            break;
                        ++z;
                    }
                    if (score >= tolerance){
						hit.position = i - window_positions[k];
						hit.score = score;
						ret[k].push_back(hit);
					}
                }
            }
        }
    }

	for (position_t i = n - q + 1; i < n; ++i) // possible hits for matrices shorter than q near the end of sequence
	{
		code = (code << BITSHIFT) & BITAND; // dummy character to the end of code

		if (!output[code].empty())
        {

			for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
        	{
            	if (y->full && m[y->matrix] < n - i + 1) // only sufficiently short hits are considered
            	{
                	hit.position = i;
                	hit.score = y->score;
                	ret[y->matrix].push_back(hit);
            	}
        	}
		}
	}

	return ret;
}

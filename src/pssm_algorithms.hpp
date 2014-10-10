// pssm_algorithms - collection of algorithms for finding PSSM matches from sequences
// Copyright (C) 2007-2009 Pasi Rastas, Janne Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.

#ifndef _PSSM_ALG_HPP
#define _PSSM_ALG_HPP

#include <vector>
#include <iostream>
#include <cfloat>
#include <stdint.h>
#include "seq_buffer.h"


// for data about position in a sequence
typedef long position_t;


// for bit parallel magic
typedef uint_fast32_t bits_t;

typedef std::vector<bits_t> bitArray;
typedef std::vector<bitArray> bitMatrix;


// matrix scores
typedef double score_t;
const score_t SCORE_MIN = DBL_MIN;
const score_t PVAL_DP_MULTIPLIER = 1000.0;


typedef std::vector<score_t> scoreArray;
typedef std::vector<scoreArray> scoreMatrix;


// Struct for storing data about matrix matches
struct matchData
{
    position_t position;
    score_t score;
};

typedef std::vector<int> intArray;
typedef std::vector<unsigned int> uintArray;
typedef std::vector<short> shortArray;
typedef std::vector<float> floatArray;
typedef std::vector<double> doubleArray;
typedef std::vector<intArray> intMatrix;
typedef std::vector<unsigned char> charArray;
typedef std::vector<matchData> matchArray;

doubleArray flatBG(const int size);
doubleArray bgFromSequence(SeqIterator &seq_it, const int size, const double pseudocount);
doubleArray bgFromSequence(const charArray &seq, const int size, const double pseudocount);
charArray readString(std::basic_istream<char> &in);
charArray readDNA(std::basic_istream<char> &in);
charArray readAZ(std::basic_istream<char> &in);
scoreMatrix readMatrix(std::basic_istream<char> &in);
scoreMatrix reverseComplement(scoreMatrix &mat);
scoreMatrix counts2LogOdds(const scoreMatrix &mat, const doubleArray &bg, const double ps);
score_t tresholdFromP(const scoreMatrix &mat, const doubleArray &bg, const double &p);
score_t maxScore(const scoreMatrix &mat);
doubleArray expectedDifferences(const scoreMatrix &mat, const doubleArray &bg);
void printScoreMatrix(const scoreMatrix &mat);
void printSequence(const charArray &seq, const int n);
void printMatchArray(const matchArray &match);


//seq iterator version of the algorithms. Supports sequences that doesn't fit to main memory
matchArray naiveAlgorithm(SeqIterator &seq_it, const scoreMatrix &p, const score_t tol);
matchArray permutatedLookAhead(SeqIterator &, const scoreMatrix &p, const doubleArray &bg, const score_t tol);

matchArray naiveSuperalphabetAlgorithmDNA(const int q, SeqIterator &, const scoreMatrix &p, const score_t tol);
matchArray lookaheadFiltrationDNA(const int q, SeqIterator &, const scoreMatrix &p, const doubleArray &bg, const score_t tol);
std::vector<matchArray> multipleMatrixLookaheadFiltrationDNA(const int q, SeqIterator &, const std::vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol);
std::vector<matchArray> searchDNA(SeqIterator &, const std::vector<scoreMatrix> &matrices, const scoreArray &tol, const doubleArray &bg, const int q);

//char array version of the algorithms
matchArray naiveAlgorithm(const charArray &s, const scoreMatrix &p, const score_t tol);
matchArray permutatedLookAhead(const charArray &s, const scoreMatrix &p, const doubleArray &bg, const score_t tol);
matchArray naiveSuperalphabetAlgorithm(const int q, const charArray &s, const scoreMatrix &p, const score_t tol);
matchArray naiveSuperalphabetAlgorithmDNA(const int q, const charArray &s, const scoreMatrix &p, const score_t tol);
matchArray lookaheadFiltration(const int q, const charArray &s, const scoreMatrix &p, const doubleArray &bg, const score_t tol);
matchArray lookaheadFiltrationDNA(const int q, const charArray &s, const scoreMatrix &p, const doubleArray &bg, const score_t tol);
std::vector<matchArray> multipleMatrixLookaheadFiltration(const int q, const charArray &s, const std::vector<intMatrix> &matrices, const doubleArray &bg, const intArray &tol);
std::vector<matchArray> multipleMatrixLookaheadFiltrationDNA(const int q, const charArray &s, const std::vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol);


#endif

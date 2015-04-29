#ifndef MOODS_H
#define MOODS_H


#include <vector>
#include <cmath>
#include <cstdint>
#include <limits>


// for data about position in a sequence
// TODO: replace with size_t 
typedef unsigned long long position_t;

// for bit parallel magic
typedef unsigned long bits_t;

// matrix scores
typedef std::vector<std::vector<double> > score_matrix;

#endif
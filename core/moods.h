#ifndef MOODS_H
#define MOODS_H


#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <limits>

// for bit parallel magic
typedef uint_fast32_t bits_t;

// matrix scores
typedef std::vector<std::vector<double> > score_matrix;


// structs for passing complex information in/out

// motif matches
struct match
{
    size_t pos;
    double score;
};

// sequence variants, e.g. SNPs and indels
struct variant
{
    size_t start_pos;
    size_t end_pos;
    std::string modified_seq;
};

// motif matches with variant information
struct match_with_variant
{
    size_t pos;
    long variant_of_pos;
    std::vector<variant> variants;
    double score;
};


#endif
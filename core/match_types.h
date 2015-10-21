#ifndef MATCH_TYPES_H
#define MATCH_TYPES_H

#include <vector>
#include <string>

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
    double score;
    std::vector<size_t> variants;
    bool starts_in_seq;
};


#endif
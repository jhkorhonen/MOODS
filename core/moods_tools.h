#ifndef MOODS_TOOLS_H
#define MOODS_TOOLS_H

#include "moods.h"
#include "match_types.h"

namespace MOODS { namespace tools{
    // background functions
    std::vector<double> flat_bg(const unsigned int alphabet_size);
    std::vector<double> bg_from_sequence_dna(const std::string &seq, const double ps);
    
    std::vector<MOODS::variant> snp_variants(const std::string &seq);
    
    // matrix transformations
    score_matrix reverse_complement(const score_matrix &mat);
    score_matrix log_odds(const score_matrix &mat, const std::vector<double> &bg, const double ps);
    score_matrix log_odds(const score_matrix &mat, const std::vector<double> &bg, const double ps, const double log_base);
    
    // threshold from p
    double threshold_from_p(const score_matrix &mat, const std::vector<double> &bg, const double &p, const double &precision);
    
    // min / max
    double max_score(const score_matrix &mat);
    double min_score(const score_matrix &mat);
    double min_delta(const score_matrix &mat);

    // high-order versions
    double max_score(const score_matrix &mat, const size_t a);
    double min_score(const score_matrix &mat, const size_t a);
    double threshold_from_p(const score_matrix &mat, const std::vector<double> &bg, const double &p, size_t a, const double &precision);
    score_matrix reverse_complement(const std::vector<std::vector<double>> &mat, size_t a);
    score_matrix log_odds(const score_matrix &mat, const std::vector<std::vector<double> >& low_order_terms,
                          const std::vector<double> &bg, const double ps, const size_t a);
    score_matrix log_odds(const score_matrix &mat, const std::vector<std::vector<double> >& low_order_terms,
                          const std::vector<double> &bg, const double ps, const size_t a, const double log_base);
    
}}

#endif
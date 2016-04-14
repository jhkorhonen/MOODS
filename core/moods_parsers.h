#ifndef MOODS_PARSERS_H
#define MOODS_PARSERS_H

#include "moods.h"

namespace MOODS { namespace parsers{

    score_matrix pfm(const std::string& filename);
    score_matrix pfm_to_log_odds(const std::string& filename, const std::vector<double> &bg, const double pseudocount, const double log_base = -1);
    
    score_matrix adm_1o_terms(const std::string& filename, const size_t a = 4);
    score_matrix adm_0o_terms(const std::string& filename, const size_t a = 4);
    
    score_matrix adm_to_log_odds(const std::string& filename, const std::vector<double> &bg, const double pseudocount, const size_t a = 4, const double log_base = -1);

}}
#endif
#ifndef MOODS_SCAN_H
#define MOODS_SCAN_H

#include <string>

#include "moods.h"

namespace MOODS { namespace scan{
    
    struct match
    {
        position_t pos;
        double score;
    };
    
    std::vector< std::vector<match> > scan_dna(const std::string& seq, const std::vector<score_matrix>& matrices, const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size );
    
}}



#endif
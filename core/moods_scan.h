#ifndef MOODS_SCAN_H
#define MOODS_SCAN_H

#include "moods.h"

namespace MOODS { namespace scan{
    
    struct scanner_output
    {
        double score;
        unsigned int matrix;
        bool full;
    };
    
    struct match
    {
        position_t pos;
        double score;
    };
    

    class Scanner {
        std::vector<Motif> motifs;
        std::vector<std::vector<scanner_output> > window_hits;
        unsigned int a;
        unsigned int av;
        unsigned int l;
    public:
        Scanner(const std::vector<score_matrix>& matrices, const std::vector<double> thresholds,  const std::vector<double> bg);
        std::vector<std::vector<match> > scan(std::vector<unsigned char>& seq);
    };
    
    
}}

#endif
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
    

    class Scanner {
        std::vector<Motif> motifs;
        std::vector<std::vector<scanner_output> > window_hits;
        unsigned int a;
        unsigned int av;
        unsigned int l;
    public:
        Scanner(const std::vector<score_matrix>& matrices, const std::vector<double> thresholds,  const std::vector<double> bg);
    };
    
    
}}

#endif
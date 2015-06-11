#ifndef MOODS_SCANNER_H
#define MOODS_SCANNER_H

#include "moods.h"
#include "motif.h"
#include "moods_scan.h"
#include "moods_misc.h"

namespace MOODS { namespace scan{
    
    struct scanner_output
    {
        double score;
        std::size_t matrix;
        bool full;
    };
    
    
    class Scanner {
        std::vector<Motif> motifs;
        std::vector<std::vector<scanner_output> > window_hits;
        unsigned int a;
        unsigned int l;
    public:
        Scanner(const std::vector<Motif>& matrices, unsigned int alphabet_size, unsigned int window_size);
        std::vector<std::vector<match> > scan(const misc::seq_internal& s);
    };
    
    
}}

#endif
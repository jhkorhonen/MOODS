#ifndef MOODS_SCANNER_H
#define MOODS_SCANNER_H

#include "moods.h"
#include "motif.h"
#include "moods_scan.h"
#include "moods_misc.h"

#include <memory>

namespace MOODS { namespace scan{
    
    struct scanner_output
    {
        double score;
        std::size_t matrix;
        bool full;
    };
    
    class Scanner {
    public:
        Scanner(unsigned int window_size);
        Scanner(unsigned int window_size, const std::vector<std::string>& alphabet);

        // void set_motifs(const std::vector<MOODS::scan::Motif>& motifs);
        void set_motifs(const std::vector<score_matrix>& matrices,
                        const std::vector<double>& bg,
                        const std::vector<double> thresholds);

        std::vector<std::vector<scan::match> > scan(const std::string& s);

        // std::vector<std::vector<scan::match> > scan(const std::string& s, size_t max_hits);

    private:
        // std::vector<MOODS::scan::Motif> motifs;
        std::vector<std::unique_ptr<MOODS::scan::Motif>> motifs;
        std::vector<std::vector<scanner_output>> window_hits;
        unsigned int a;
        unsigned int l;
        std::vector<unsigned char> alphabet_map;
        bool initialised = false;

        void initialise_hit_table();
        std::vector<size_t> preprocess_seq(const std::string& s); 
        template<typename T> void process_matches(const std::string& s, const T& match_handler);        
    };
}}

#endif
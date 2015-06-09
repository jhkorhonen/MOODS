

#include "moods.h"
#include "moods_scan.h"
#include "moods_misc.h"
#include "motif.h"
#include "scanner.h"
#include <iostream>

using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{
    
    vector< vector< scan::match> > scan_dna(const string& seq, const vector<score_matrix>& matrices, const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    {
        
        vector<Motif> motifs;
        
        std::cout << "Building motif objects\n";
        
        for (size_t i = 0; i < matrices.size(); ++i){
            motifs.emplace_back( matrices[i], bg, window_size, thresholds[i] );
        }
        
        std::cout << "Building scanner object\n";
        
        Scanner scanner(motifs, 4, window_size);
        
        std::cout << "Parsing sequence\n";
        
        misc::seq_internal s = misc::string_to_seq_dna(seq);
        
        std::cout << "Scanning\n";
        
        return scanner.scan(s);
    }
    
}}




#include "moods.h"
#include "moods_scan.h"
#include "moods_misc.h"
#include "motif.h"
#include "scanner.h"
#include <iostream>
#include <iomanip>


using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{
    
    vector< vector< scan::match> > scan_dna(const string& seq, const vector<score_matrix>& matrices, const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    // unsigned int scan_dna(const string& seq, const vector<score_matrix>& matrices, const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    {
        
        clock_t total = clock();

        
        vector<Motif> motifs;
        
        std::cerr << "Building motif objects\n";
        
        clock_t start = clock();
        
        for (size_t i = 0; i < matrices.size(); ++i){
            motifs.emplace_back( matrices[i], bg, window_size, thresholds[i] );
        }
        
        std::cerr << "Building scanner object\n";
        
        Scanner scanner(motifs, window_size);
        
        std::cerr << "Preprocessing motifs: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        std::cerr << "Parsing sequence\n";
        start = clock();
        
        // misc::seq_internal s = misc::string_to_seq_dna(seq);
        //
        // std::cerr << "Preprocessing sequence: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";

        std::cerr << "Scanning\n";        
        start = clock();
        
        // for (size_t i = 0; i < s.starts.size(); ++i){
        //     std::cerr << s.starts[i] << "\n";
        // }
        //
        // for (size_t i = 0; i < s.ends.size(); ++i){
        //     std::cerr << s.ends[i] << "\n";
        // }
        //
        
        auto results = scanner.scan(seq);
        
        std::cerr << "Scanning: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        std::cerr << "Total: " << std::setprecision(5) << (double)(clock() - total)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        
        return results;
    }

    std::vector< std::vector< match> > scan(const std::string& seq, const std::vector<score_matrix>& matrices, const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size, const std::vector<std::string>& alphabet)
    {
        clock_t total = clock();
                
        vector<Motif> motifs;
        
        std::cerr << "Building motif objects\n";
        
        clock_t start = clock();
        
        for (size_t i = 0; i < matrices.size(); ++i){
            motifs.emplace_back( matrices[i], bg, window_size, thresholds[i] );
        }
        
        std::cerr << "Building scanner object\n";
        
        Scanner scanner(motifs, window_size, alphabet);
        
        std::cerr << "Preprocessing motifs: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        std::cerr << "Parsing sequence\n";
        start = clock();
        
        // misc::seq_internal s = misc::string_to_seq_dna(seq);
        //
        // std::cerr << "Preprocessing sequence: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";

        std::cerr << "Scanning\n";        
        start = clock();
        
        // for (size_t i = 0; i < s.starts.size(); ++i){
        //     std::cerr << s.starts[i] << "\n";
        // }
        //
        // for (size_t i = 0; i < s.ends.size(); ++i){
        //     std::cerr << s.ends[i] << "\n";
        // }
        //
        
        auto results = scanner.scan(seq);
        
        std::cerr << "Scanning: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        std::cerr << "Total: " << std::setprecision(5) << (double)(clock() - total)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        
        return results;
    }
    
}}


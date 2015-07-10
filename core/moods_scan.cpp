

#include "moods.h"
#include "moods_scan.h"
#include "moods_misc.h"
#include "motif.h"
#include "scanner.h"
// #include <iostream>
// #include <iomanip>


using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{

    vector< vector< scan::match> > scan_dna(const string& seq, const vector<score_matrix>& matrices, const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    {
        
        // clock_t total = clock();
        // clock_t start = clock();
        
        // std::cerr << "Building scanner object\n";
        
        Scanner scanner(window_size);

        scanner.set_motifs(matrices, bg, thresholds);
        
        // std::cerr << "Preprocessing motifs: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        // std::cerr << "Scanning\n";        
        // start = clock();

        auto results = scanner.scan(seq);
        
        // std::cerr << "Scanning: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        // std::cerr << "Total: " << std::setprecision(5) << (double)(clock() - total)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        
        return results;
    }

    std::vector< std::vector< match> > scan(const std::string& seq, const std::vector<score_matrix>& matrices, const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size, const std::vector<std::string>& alphabet)
    {
        // clock_t total = clock();
        // clock_t start = clock();
        
        // std::cerr << "Building scanner object\n";
        
        Scanner scanner(window_size, alphabet);

        scanner.set_motifs(matrices, bg, thresholds);
        
        // std::cerr << "Preprocessing motifs: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";

        // std::cerr << "Scanning\n";        
        // start = clock();
        
        auto results = scanner.scan(seq);
        
        // std::cerr << "Scanning: " << std::setprecision(5) << (double)(clock() - start)/((double)CLOCKS_PER_SEC) <<     "\n";
        // std::cerr << "Total: " << std::setprecision(5) << (double)(clock() - total)/((double)CLOCKS_PER_SEC) <<     "\n";
        
        return results;
    }
    
}}


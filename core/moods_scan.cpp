

#include "moods.h"
#include "moods_scan.h"
#include "moods_misc.h"
#include "motif.h"
#include "scanner.h"

using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{

    vector< vector<match> > scan_dna(const string& seq, const vector<score_matrix>& matrices,
                                            const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    {
        
        Scanner scanner(window_size);
        scanner.set_motifs(matrices, bg, thresholds);

        auto results = scanner.scan(seq);        
        return results;
    }

    std::vector< std::vector< match> > scan(const std::string& seq, const std::vector<score_matrix>& matrices,
                                            const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size,
                                            const std::vector<std::string>& alphabet)
    {

        Scanner scanner(window_size, alphabet);
        scanner.set_motifs(matrices, bg, thresholds);
        
        auto results = scanner.scan(seq);
        return results;
    }

    // void variant_matches_recursive(std::vector<std::vector<match>> results, ...,
    //                                scanner& scanner, const std::string& seq, const std::vector<variant> variants){
        
    // }

    std::vector<std::vector<match_with_variant>> variant_matches(MOODS::scan::Scanner& scanner, const std::string& seq, const std::vector<variant> variants){
        return vector<vector<match_with_variant>>();
    }
    
}}




#include "moods.h"
#include "moods_scan.h"

using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{
    
    vector< vector<match> > scan_dna(const string& seq, const vector<score_matrix>& matrices, const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    {
        // ...
    }
    
}}


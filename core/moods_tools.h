
#include <vector>
#include <iostream>
#include <cfloat>
#include <stdint.h>

#include "moods.h"

namespace MOODS { namespace tools{


// background functions
std::vector<double> flat_bg(const unsigned int alphabet_size);
std::vector<double> bg_from_sequence(const std::vector<unsigned char> &seq, const int alphabet_size, const double ps);

// matrix transformations
score_matrix reverse_complement(const score_matrix &mat);
score_matrix log_odds(const score_matrix &mat, const std::vector<double> &bg, const double ps);
score_matrix log_odds(const score_matrix &mat, const std::vector<double> &bg, const double ps, const double log_base);

// min / max
double max_score(const score_matrix &mat);
double min_score(const score_matrix &mat);


}}
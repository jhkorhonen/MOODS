
#include <vector>
#include <string>

#include "match_types.h"

namespace MOODS {
    
    variant::variant(){}
    variant::variant(size_t s, size_t e, std::string seq): start_pos(s), end_pos(e), modified_seq(seq) {}
    
}
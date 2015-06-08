#ifndef MOODS_MISC_H
#define MOODS_MISC_H

#include "moods.h"

namespace MOODS { namespace misc{
    
    struct seq_internal
    {
        std::vector<unsigned char> seq;
        std::vector<std::size_t> starts;
        std::vector<std::size_t> ends;
    };
    
    unsigned int shift(unsigned int a);
    bits_t mask(unsigned int a);
    
    seq_internal string_to_seq_dna(const std::string& s);
}}



#endif
#ifndef MOODS_MISC_H
#define MOODS_MISC_H

#include "moods.h"

namespace MOODS { namespace misc{
    
    unsigned int shift(unsigned int a);
    bits_t mask(unsigned int a);
    unsigned int q_gram_size(size_t rows, unsigned int a);

    std::vector<size_t> preprocess_seq(const std::string& s, size_t a, const std::vector<unsigned char>& alphabet_map); 

}}



#endif
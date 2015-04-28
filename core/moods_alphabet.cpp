

#include "moods.h"
#include "moods_alphabet.h"


namespace MOODS { namespace alphabet{

    // basically base-2 logarithm of a, rounded up
    unsigned int shift(unsigned int a)
    {
        unsigned int s = 0;
        unsigned int b = 1;
        
        while (b < a){
            s += 1;
            b =<< 1;
        }
        return s;
    }
        
    bits_t mask(unsigned int a){
        bits_t b = 1;
        
        while (b < a){
            b =<< 1;
        }
        return b-1;
    }

}}
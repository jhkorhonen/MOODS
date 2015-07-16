#include "moods.h"
#include "moods_misc.h"


using std::vector;
using std::size_t;

namespace MOODS { namespace misc{

    // basically base-2 logarithm of a, rounded up
    unsigned int shift(unsigned int a)
    {
        unsigned int s = 0;
        unsigned int b = 1;
        
        while (b < a){
            s += 1;
            b = b << 1;
        }
        return s;
    }
        
    bits_t mask(unsigned int a){
        bits_t b = 1;
        
        while (b < a){
            b = b << 1;
        }
        return b-1;
    }

    unsigned int q_gram_size(size_t rows, unsigned int a){
        unsigned int q = 0;
        unsigned int s = 1;
        

        while (s < rows){
            q += 1;
            s *= a;
        }
        return q;   
    }

}}
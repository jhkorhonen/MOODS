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

    // checks a sequence for non-scan regions and returns the corresponding bounds
    std::vector<size_t> preprocess_seq(const std::string& s, size_t a, const std::vector<unsigned char>& alphabet_map){
        
        vector<size_t> bounds;
        
        bool scannable = false;
        unsigned char c;
        
        for (size_t i = 0; i < s.size(); ++i){
            c = alphabet_map[(unsigned char)s[i]];
            
            if (c < a){
                if (!scannable){
                    scannable = true;
                    bounds.push_back(i);
                }
            }
            else {
                if (scannable){
                    scannable = false;
                    bounds.push_back(i);
                }
            }
        }
        if (scannable){
            bounds.push_back(s.size());
        }
        
        return bounds;
        
    }
}}


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
    
    // maps ACGT string to the interal presentation 
    seq_internal string_to_seq_dna(const std::string& s){
        
        // TODO: re-write this part
        vector<char> m(4, sizeof(char));
        
        m['a'] = 0;
        m['A'] = 0;
        
        m['c'] = 1;
        m['C'] = 1;
        
        m['g'] = 2;
        m['G'] = 2;
        
        m['t'] = 3;
        m['T'] = 3;
        
        
        seq_internal ret = {
            vector<unsigned char> (4, s.size()),
            vector<size_t>(),
            vector<size_t>()
        };
        
        
        // TODO: re-write this part also...
        bool scannable = false; 
        for (size_t i = 0; i < s.size(); ++i){
            unsigned char c = m[s[i]];
            
            if (c < 4){
                ret.seq[i] = c;
                if (!scannable){
                    scannable = true;
                    ret.starts.push_back(i);
                }
            }
            else {
                if (scannable){
                    scannable = false;
                    ret.ends.push_back(i);
                }
            }
        }
        if (scannable){
            ret.starts.push_back(s.size());
        }
        
        return ret;
    }
    

}}


#include "moods.h"
#include "moods_misc.h"
#include <iostream>


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
        vector<unsigned char> m(127, 4);
        
        m['a'] = 0;
        m['A'] = 0;
        
        m['c'] = 1;
        m['C'] = 1;
        
        m['g'] = 2;
        m['G'] = 2;
        
        m['t'] = 3;
        m['T'] = 3;        
        
        vector<unsigned char> sequence(s.size(), 4);
        vector<size_t> start_pos; 
        vector<size_t> end_pos;
        
        
        // TODO: re-write this part also...
        bool scannable = false; 
        for (size_t i = 0; i < s.size(); ++i){
            unsigned char c = m[(int)s[i]];
            
            if (c < 4){
                sequence[i] = c;
                if (!scannable){
                    scannable = true;
                    start_pos.push_back(i);
                    std::cout << "start " << i << "\n";
                }
            }
            else {
                if (scannable){
                    scannable = false;
                    end_pos.push_back(i);
                    std::cout << "end " << i << "\n";
                }
            }
        }
        if (scannable){
            end_pos.push_back(s.size());
            std::cout << "end " << s.size() << "\n";
        }
        
        seq_internal ret;
        ret.seq = sequence;
        ret.starts = start_pos;
        ret.ends = end_pos;
        
        return ret;
    }
    

}}
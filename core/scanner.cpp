// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include <utility>
#include <tuple>
#include <iostream>

#include "moods.h"
#include "scanner.h"
#include "motif.h"
#include "moods_misc.h"

using std::vector;
using std::size_t;


namespace MOODS { namespace scan{

    Scanner::Scanner(const std::vector<Motif>& matrices, unsigned int alphabet_size, unsigned int window_size)
    {
        motifs = matrices;
        a = alphabet_size;
        l = window_size;
        
        // for (unsigned int k; k < matrices.size(), ++k)
        // {
        //     Motif m = Motif(matrices[k], bg, l, thresholds[k]);
        //     motifs.push_back(m);
        // }
        
        const bits_t SHIFT = MOODS::misc::shift(a);
        const bits_t CODE_SIZE = 1 << (SHIFT * l);
        
        window_hits = vector<vector<scanner_output> >(CODE_SIZE, vector<scanner_output>());
        
        for (bits_t code = 0; code < CODE_SIZE; ++code)
        {
            for (size_t k = 0; k < motifs.size(); ++k)
            {
                std::cout << code << " " << k << "\n";
                double score;
                bool match;
                
                std::tie(match,score) = motifs[k].window_match(code, SHIFT);
                
                if (match)
                {
                    scanner_output op = {score, k, motifs[k].size() <= l};
                    window_hits[code].emplace_back(op);
                }
            }
        }
    }
    
    std::vector<std::vector<match> > Scanner::scan(const misc::seq_internal& s)
    {
        const bits_t SHIFT = MOODS::misc::shift(a);
        const bits_t MASK = (1 << (SHIFT * l)) - 1;
        
        vector<vector<match> > ret(motifs.size(), vector<match>());
        
        const vector<unsigned char>& seq = s.seq;
        
        std::cout << s.starts.size() << " " << s.ends.size() << "\n";
        
        // Scanning
        for (size_t seq_i = 0; seq_i < s.starts.size(); ++seq_i){
            size_t start = s.starts[seq_i];
            size_t end = s.ends[seq_i];
            
            
            // sequence is very short
            if (end - start < l){
                std::cout << "short\n";
                
                bits_t code = 0;
                for (size_t i = start; i < end; ++i)
                    code = (code << SHIFT) + seq[i];
                
                for (size_t i = end - start; i < l - 1; ++ i){
                    code = (code << SHIFT) & MASK;  // dummy character to the end of code
                }
                
                
                for (size_t i = start; i < end; ++i)
                    if (!window_hits[code].empty())
                    {
                        code = (code << SHIFT) & MASK;  // dummy character to the end of code
                        
                        for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                        {
                            if (y->full && motifs[y->matrix].size() <= end - i) // only sufficiently short hits are considered
                            {
                                ret[y->matrix].push_back(match{i,y->score});
                            }
                        }
                    }
                }
            // sequence is long enough that we have at least one proper scanning step
            else {
                
                std::cout << "long\n";
                // Initialise scanner state
                bits_t code = 0;
                for (size_t i = start; i < start + l - 1; ++i){
                    std::cout << code << " " << i << "\n";
                    code = (code << SHIFT) + seq[i];
                }
                
                // Actual scanning for the 'middle' of the sequence
                for (size_t i = start; i < end - l + 1; ++i)
                {
                    code = ((code << SHIFT) + seq[i + l - 1]) & MASK;
                    std::cout << code << " " << i << "\n";
                
                    if (!window_hits[code].empty())
                    {
                        for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                        {
                            if (y->full) // A Hit for a matrix of length <= q
                            {
                                ret[y->matrix].push_back(match{i,y->score});
                                continue;
                            }
                            if (i - start >= motifs[y->matrix].window_pos() && i + motifs[y->matrix].size() - motifs[y->matrix].window_pos() <= end) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                            {
                                std::cout << start << " " << i << " "  << motifs[y->matrix].window_pos() <<  "\n";
                                double score = motifs[y->matrix].check_hit(seq, i, y->score);
                                if (score >= motifs[y->matrix].threshold()){
                                    ret[y->matrix].push_back(match{i - motifs[y->matrix].window_pos(),score});
                                }
                            }
                        }
                    }
                }
                
                std::cout << "tail\n";
            
                // possible hits for matrices shorter than l near the end of current interval
                for (size_t i = end - l + 1; i < end; ++i)
                {
                    code = (code << SHIFT) & MASK;  // dummy character to the end of code
                    std::cout << code << "\n";
            
                    if (!window_hits[code].empty())
                    {
                
                        for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                        {
                            if (y->full && motifs[y->matrix].size() < end - i) // only sufficiently short hits are considered
                            {
                                ret[y->matrix].push_back(match{i,y->score});
                            }
                        }
                    }
                }
            } 
        }
        
        return ret;
    }


} // namespace scan
} // namespace MOODS
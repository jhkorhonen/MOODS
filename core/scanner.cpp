// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include <utility>
#include <tuple>

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
                double score;
                bool match;
                
                std::tie(match,score) = motifs[k].window_match(code, SHIFT);
                
                if (match)
                {
                    scanner_output op = {score, k, motifs[k].size() <= l};
                    window_hits[code].push_back(op);
                }
            }
        }
    }
    
    std::vector<std::vector<match> > Scanner::scan(vector<unsigned char>& seq)
    {
        const bits_t SHIFT = MOODS::misc::shift(a);
        const bits_t MASK = (1 << (SHIFT * l)) - 1;
        
        // Initialise scanner state
        bits_t code = 0;
        for (size_t i = 0; i < l - 1; ++i)
            code = (code << SHIFT) + seq[i];
        
        
        vector<vector<match> > ret(vector<match>(), motifs.size());
        
        
        // Scanning
        for (size_t i = 0; i < n - l + 1; ++i)
        {
        	code = ((code << SHIFT) + seq[i + l - 1]) & MASK;

            if (!window_hits[code].empty())
            {
                for (vector<OutputListElementMulti>::iterator y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                {
                    if (y->full) // A Hit for a matrix of length <= q
                    {
                        ret[y->matrix].emplace_back(i,y->score);
                        continue;
                    }
                    if (i >= motifs[y->matrix].window_pos() && i + motifs[y->matrix].size() - motifs[y->matrix].window_pos() <= n) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                    {
                        double score = motifs[y->matrix].check_hit(seq, i, y->score);
                        if (score >= motifs[y->matrix].threshold()){
                            ret[y->matrix].emplace_back(i,score);
                        }
                    }
                }
            }
        }
        
        for (size_t i = n - l + 1; i < n; ++i) // possible hits for matrices shorter than l near the end of sequence
        {
        	code = (code << SHIFT) & MASK;  // dummy character to the end of code
            
            if (!window_hits[code].empty())
            {
                
                for (vector<OutputListElementMulti>::iterator y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                {
                    if (y->full && motif[y->matrix].size() < n - i + 1) // only sufficiently short hits are considered
                    {
                        ret[y->matrix].emplace_back(i,y->score);
                    }
                }
            }
        }
        
        return ret;
    }


} // namespace scan
} // namespace MOODS
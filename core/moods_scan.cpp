// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include "moods.h"
#include "moods_scan.h"
#include "moods_alphabet.h"

using std::vector;


namespace MOODS { namespace scan{

    Scanner::Scanner(const std::vector<score_matrix>& matrices, const std::vector<double> thresholds,  const std::vector<double> bg, unsigned int window_size)
    {
        motifs = vector<Motif>();
        a = bg.size();        
        l = window_size;
        
        for (unsigned int k; k < matrices.size(), ++k)
        {
            Motif m = Motif(matrices[k], bg, l, thresholds[k]);
            motifs.push_back(m);
        }
        
        const bits_t SHIFT = MOODS::alphabet::shift(a);
        const bits_t CODE_SIZE = 1 << (SHIFT * l);
        
        window_hits = vector<vector<scanner_output> >(CODE_SIZE, vector<scanner_output>());
        
        for (bits_t code = 0; code < CODE_SIZE; ++code)
        {
            for (unsigned int k = 0; k < motifs.size(); ++k)
            {
                if (motifs[k].window_match(code, SHIFT))
                {
                    scanner_output op = {motifs[k].window_score(code, SHIFT), k, motifs[k].size() <= l};
                    window_hits[code].push_back(op);
                }
            }
        }
    }


} // namespace scan
} // namespace MOODS
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
#include "scanner_tools.h"

using std::vector;
using std::size_t;


namespace MOODS { namespace scan{

	MOODS::scan::Scanner build_scanner(const std::vector<score_matrix>& matrices, const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size)
    {
        vector<Motif> motifs;
        
        for (size_t i = 0; i < matrices.size(); ++i){
            motifs.emplace_back( matrices[i], bg, window_size, thresholds[i] );
        }

        Scanner scanner(motifs, window_size);

        return scanner;
    }

} // namespace scan
} // namespace MOODS

#ifndef MOODS_SCANNER_TOOLS_H
#define MOODS_SCANNER_TOOLS_H

#include "moods.h"
#include "motif.h"
#include "scanner.h"
#include "moods_scan.h"
#include "moods_misc.h"


namespace MOODS { namespace scan{

    MOODS::scan::Scanner build_scanner(const std::vector<score_matrix>& matrices,
                      const std::vector<double>& bg,
                      const std::vector<double> thresholds,
                      unsigned int window_size);

    
    
}}

#endif
#ifndef MOODS_PARSERS_H
#define MOODS_PARSERS_H

#include "moods.h"

namespace MOODS { namespace parsers{

    std::vector<std::vector<double> > pfm(const std::string& filename);
    std::vector<std::vector<double> > pfm_log_odds(const std::string& filename, const std::vector<double> &bg, double pseudocount);
    std::vector<std::vector<double> > pfm_log_odds_rc(const std::string& filename, const std::vector<double> &bg, double pseudocount);


    std::vector<std::vector<double> > adm_1o_terms(const std::string& filename, size_t a = 4);
    std::vector<std::vector<double> > adm_0o_terms(const std::string& filename, size_t a = 4);

    std::vector<std::vector<double> > adm_log_odds(const std::string& filename, const std::vector<double> &bg, double pseudocount, size_t a = 4);
    std::vector<std::vector<double> > adm_log_odds_rc(const std::string& filename, const std::vector<double> &bg, double pseudocount, size_t a = 4);

}}



#endif
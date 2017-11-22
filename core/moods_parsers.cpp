// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include "moods.h"
#include "moods_parsers.h"
#include "moods_tools.h"
#include "moods_misc.h"

#include <climits>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

using std::vector;
using std::size_t;
using std::string;

namespace MOODS { namespace parsers{

    vector<vector<double>> read_table(const string& filename){
        std::ifstream ifs (filename, std::ifstream::in);
        vector<vector<double>> mat;     

        string line;
        while (std::getline(ifs, line)){
            vector<double> row;
            std::istringstream iss(line);
            std::copy(std::istream_iterator<double>(iss),
                      std::istream_iterator<double>(),
                      std::back_inserter(row));
            
            // let's not put any empty lines in there 
            if (row.size() >= 1){
                mat.push_back(row);
            }
        }

        return mat;

    }

    score_matrix pfm(const string& filename){
        score_matrix mat = read_table(filename);
        
        size_t a = mat.size();
        size_t n = mat[0].size();
        
        if (a == 0 or n == 0){
            return score_matrix();
        }
        
        for (size_t i = 0; i < a; ++i){
            if (mat[i].size() != n){
                return score_matrix();
            }
        }
        
        return mat;
    }

    score_matrix pfm_to_log_odds(const string& filename, const vector<double> &bg, const double pseudocount, const double log_base){
        score_matrix mat = read_table(filename);
        
        size_t a = mat.size();
        size_t n = mat[0].size();
        
        if (a == 0 or n == 0){
            return score_matrix();
        }
        
        for (size_t i = 0; i < a; ++i){
            if (mat[i].size() != n){
                return score_matrix();
            }
        }
        
        if (log_base < 0){
            return tools::log_odds(mat, bg, pseudocount);
        }
        else {
            return tools::log_odds(mat, bg, pseudocount, log_base);
        }
    }
    
    
    vector<vector<double>> read_and_check_adm(const string& filename, const size_t a){
        vector<vector<double>> adm = read_table(filename);
        
        if (adm.size() != a * a + a){
            return score_matrix();
        }
        
        size_t n = adm[0].size();
        
        for (size_t i = 0; i < a * a; ++i){
            if (adm[i].size() != n){
                return score_matrix();
            }
        }
        
        for (size_t i = a * a + 1; i < a * a + a; ++i){
            if (adm[i].size() == 0){ // we only need 0-order term for the first row
                return score_matrix();
            }
        }
        
        return adm;
    }

    vector<vector<double>> adm_1o_terms(const string& filename, const size_t a){
        vector<vector<double>> adm = read_and_check_adm(filename,a);
        if (adm.size() == 0){
            return adm;
        }
        
        vector<vector<double>> ret;
        
        for (size_t i = 0; i < a * a; ++i){
            ret.push_back(adm[i]);
        }

        return ret;
    }

    vector<vector<double>> adm_0o_terms(const string& filename, const size_t a){
        vector<vector<double>> adm = read_and_check_adm(filename,a);
        
        if (adm.size() == 0){
            return adm;
        }
        
        vector<vector<double>> ret;
        for (size_t i = a * a; i < a * a + a; ++i){
            ret.push_back(adm[i]);
        }

        return ret;
    }

    score_matrix adm_to_log_odds(const string& filename, const vector<double> &bg,
                                        const double pseudocount, const size_t a, const double log_base){
        vector<vector<double>> adm = read_and_check_adm(filename,a);
        
        if (adm.size() == 0){
            return adm;
        }
        
        vector<vector<double>> mat;
        for (size_t i = 0; i < a * a; ++i){
            mat.push_back(adm[i]);
        }
        vector<vector<double>> zero_terms (1, vector<double>(a, 0));
        for (size_t i = 0; i < a; ++i){
            zero_terms[0][i] = adm[a * a + i][0];
        }
        return tools::log_odds(mat, zero_terms, bg, pseudocount, a);
    }

} // namespace tools
} // namespace MOODS
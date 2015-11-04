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

            mat.push_back(row);
        }

        return mat;

    }

    vector<vector<double>> pfm(const string& filename){
        // pfm files don't have any metadata so...
        return read_table(filename);
    }

    vector<vector<double>> pfm_log_odds(const string& filename, const vector<double> &bg, double pseudocount){
        // pfm files don't have any metadata so...
        return tools::log_odds(read_table(filename), bg, pseudocount);
    }

    vector<vector<double>> pfm_log_odds_rc(const string& filename, const vector<double> &bg, double pseudocount){
        // pfm files don't have any metadata so...
        return tools::log_odds(tools::reverse_complement(read_table(filename)), bg, pseudocount);
    }



    vector<vector<double>> adm_1o_terms(const string& filename, size_t a){
        vector<vector<double>> adm = read_table(filename);
        vector<vector<double>> ret;
        for (size_t i = 0; i < a * a; ++i){
            ret.push_back(adm[i]);
        }

        return ret;
    }


    vector<vector<double>> adm_0o_terms(const string& filename, size_t a){
        vector<vector<double>> adm = read_table(filename);
        vector<vector<double>> ret;
        for (size_t i = a * a; i < a * a + a; ++i){
            ret.push_back(adm[i]);
        }

        return ret;
    }

    vector<vector<double>> adm_log_odds(const string& filename, const vector<double> &bg,
                                        double pseudocount, size_t a){
        vector<vector<double>> adm = read_table(filename);
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

    vector<vector<double>> adm_log_odds_rc(const string& filename, const vector<double> &bg,
                                        double pseudocount, size_t a){
        vector<vector<double>> adm = read_table(filename);
        vector<vector<double>> mat;
        for (size_t i = 0; i < a * a; ++i){
            mat.push_back(adm[i]);
        }
        vector<vector<double>> zero_terms (1, vector<double>(a, 0));
        for (size_t i = 0; i < a; ++i){
            zero_terms[0][i] = adm[a * a + i][0];
        }
        return tools::log_odds_rc(mat, zero_terms, bg, pseudocount, a);
    }

} // namespace tools
} // namespace MOODS
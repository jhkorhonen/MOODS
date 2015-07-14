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
using std::unique_ptr;


namespace MOODS { namespace scan{

    Scanner::Scanner(unsigned int window_size){

        a = 4;
        l = window_size;

        alphabet_map = vector<unsigned char>(256, 4);
        
        alphabet_map[(unsigned char)'a'] = 0;
        alphabet_map[(unsigned char)'A'] = 0;
        
        alphabet_map[(unsigned char)'c'] = 1;
        alphabet_map[(unsigned char)'C'] = 1;
        
        alphabet_map[(unsigned char)'g'] = 2;
        alphabet_map[(unsigned char)'G'] = 2;
        
        alphabet_map[(unsigned char)'t'] = 3;
        alphabet_map[(unsigned char)'T'] = 3; 
    }

    Scanner::Scanner(unsigned int window_size, const std::vector<std::string>& alphabet)
    {
        a = alphabet.size();
        l = window_size;

        alphabet_map = vector<unsigned char>(256, a);

        for (size_t i = 0; i < alphabet.size(); ++i){
            for (size_t j = 0; j < alphabet[i].size(); ++j){
                alphabet_map[(unsigned int)alphabet[i][j]] = i;
            }
        }

        this->initialise_hit_table();
    }

    // void Scanner::set_motifs(const std::vector<MOODS::scan::Motif>& motifs){
    //     this->motifs = motifs;
    //     this->initialise_hit_table();
    // }
    
    void Scanner::set_motifs(const std::vector<score_matrix>& matrices,
                        const std::vector<double>& bg,
                        const std::vector<double> thresholds){

        this->motifs = vector<unique_ptr<Motif>>();
        
        for (size_t i = 0; i < matrices.size(); ++i){
            motifs.emplace_back(new Motif0(matrices[i], bg, l, thresholds[i]));
        }
        
        this->initialise_hit_table();

    }

    void Scanner::initialise_hit_table(){

        const bits_t SHIFT = MOODS::misc::shift(a);
        const bits_t CODE_SIZE = 1 << (SHIFT * l);
        
        window_hits = vector<vector<scanner_output> >(CODE_SIZE, vector<scanner_output>());
        
        for (bits_t code = 0; code < CODE_SIZE; ++code)
        {
            for (size_t k = 0; k < motifs.size(); ++k)
            {
                double score;
                bool match;
                
                std::tie(match,score) = motifs[k]->window_match(code, SHIFT);
                
                if (match)
                {
                    scanner_output op = {score, k, motifs[k]->size() <= l};
                    window_hits[code].emplace_back(op);
                }
            }
        }

        initialised = true;
    }

    
    // checks a sequence for non-scan regions and returns the corresponding bounds
    std::vector<size_t> Scanner::preprocess_seq(const std::string& s){
        
        vector<size_t> bounds;
        
        bool scannable = false;
        unsigned char c;
        
        for (size_t i = 0; i < s.size(); ++i){
            c = alphabet_map[(unsigned char)s[i]];
            
            if (c < a){
                if (!scannable){
                    scannable = true;
                    bounds.push_back(i);
                }
            }
            else {
                if (scannable){
                    scannable = false;
                    bounds.push_back(i);
                }
            }
        }
        if (scannable){
            bounds.push_back(s.size());
        }
        
        return bounds;
        
    }
    
    std::vector<std::vector<match> > Scanner::scan(const std::string& s)
    {

        if (!initialised){
            return vector<vector<match>>(0, vector<match>());
        }

        const bits_t SHIFT = MOODS::misc::shift(a);
        const bits_t MASK = (1 << (SHIFT * l)) - 1;
        
        vector<vector<match> > ret(motifs.size(), vector<match>());
        
        vector<size_t> bounds = preprocess_seq(s);
        
        // Scanning
        for (size_t seq_i = 0; seq_i < bounds.size(); ){
            size_t start = bounds[seq_i];
            ++seq_i;
            size_t end = bounds[seq_i];
            ++seq_i;
            
            
            // sequence is very short
            if (end - start < l){
                
                bits_t code = 0;
                for (size_t i = start; i < end; ++i)
                    code = (code << SHIFT) + alphabet_map[s[i]];
                
                for (size_t i = end - start; i < l - 1; ++ i){
                    code = (code << SHIFT) & MASK;  // dummy character to the end of code
                }
                
                
                for (size_t i = start; i < end; ++i)
                    if (!window_hits[code].empty())
                    {
                        code = (code << SHIFT) & MASK;  // dummy character to the end of code
                        
                        for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                        {
                            if (y->full && motifs[y->matrix]->size() <= end - i) // only sufficiently short hits are considered
                            {
                                ret[y->matrix].push_back(match{i,y->score});
                            }
                        }
                    }
                }
            // sequence is long enough that we have at least one "proper" scanning step
            else {
                // Initialise scanner state
                bits_t code = 0;
                for (size_t i = start; i < start + l - 1; ++i){
                    code = (code << SHIFT) + alphabet_map[s[i]];
                }
                
                // Actual scanning for the 'middle' of the sequence
                for (size_t i = start; i < end - l + 1; ++i)
                {
                    code = ((code << SHIFT) + alphabet_map[s[i + l - 1]]) & MASK;

                    if (!window_hits[code].empty())
                    {
                        for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                        {
                            if (y->full) // A Hit for a matrix of length <= q
                            {
                                ret[y->matrix].push_back(match{i,y->score});
                                continue;
                            }
                            if (i - start >= motifs[y->matrix]->window_pos() && i + motifs[y->matrix]->size() - motifs[y->matrix]->window_pos() <= end) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                            {
                                double score;
                                bool m;
                                std::tie(m,score) = motifs[y->matrix]->check_hit(s, alphabet_map, i, y->score);
                                if (m){
                                    ret[y->matrix].push_back(match{i - motifs[y->matrix]->window_pos(),score});
                                }
                            }
                        }
                    }
                }
            
                // possible hits for matrices shorter than l near the end of current interval
                for (size_t i = end - l + 1; i < end; ++i)
                {
                    code = (code << SHIFT) & MASK;  // dummy character to the end of code
            
                    if (!window_hits[code].empty())
                    {
                
                        for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
                        {
                            if (y->full && motifs[y->matrix]->size() < end - i) // only sufficiently short hits are considered
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


    // -------------
    // TODO: re-write this without code duplication
    // -------------
    // std::vector<std::vector<match> > Scanner::scan(const std::string& s, size_t max_hits)
    // {

    //     if (!initialised){
    //         return vector<vector<match>>(0, vector<match>());
    //     }

    //     const bits_t SHIFT = MOODS::misc::shift(a);
    //     const bits_t MASK = (1 << (SHIFT * l)) - 1;
        
    //     vector<vector<match> > ret(motifs.size(), vector<match>());

    //     vector<size_t> hits(motifs.size(), 0);
    //     size_t matrices_left = motifs.size();
        
    //     vector<size_t> bounds = preprocess_seq(s);
        
    //     // Scanning
    //     for (size_t seq_i = 0; seq_i < bounds.size(); ){
    //         size_t start = bounds[seq_i];
    //         ++seq_i;
    //         size_t end = bounds[seq_i];
    //         ++seq_i;
            
            
    //         // sequence is very short
    //         if (end - start < l){
                
    //             bits_t code = 0;
    //             for (size_t i = start; i < end; ++i)
    //                 code = (code << SHIFT) + alphabet_map[s[i]];
                
    //             for (size_t i = end - start; i < l - 1; ++ i){
    //                 code = (code << SHIFT) & MASK;  // dummy character to the end of code
    //             }
                
                
    //             for (size_t i = start; i < end; ++i)
    //                 if (!window_hits[code].empty())
    //                 {
    //                     code = (code << SHIFT) & MASK;  // dummy character to the end of code
                        
    //                     for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
    //                     {
    //                         if (y->full && motifs[y->matrix].size() <= end - i) // only sufficiently short hits are considered
    //                         {
    //                             if (hits[y->matrix] < max_hits){
    //                                 ret[y->matrix].push_back(match{i,y->score});
    //                                 hits[y->matrix]++;
    //                                 if (hits[y->matrix] == max_hits){
    //                                     matrices_left--;
    //                                     if (matrices_left == 0){
    //                                         return ret;
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         // sequence is long enough that we have at least one "proper" scanning step
    //         else {
    //             // Initialise scanner state
    //             bits_t code = 0;
    //             for (size_t i = start; i < start + l - 1; ++i){
    //                 code = (code << SHIFT) + alphabet_map[s[i]];
    //             }
                
    //             // Actual scanning for the 'middle' of the sequence
    //             for (size_t i = start; i < end - l + 1; ++i)
    //             {
    //                 code = ((code << SHIFT) + alphabet_map[s[i + l - 1]]) & MASK;

    //                 if (!window_hits[code].empty())
    //                 {
    //                     for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
    //                     {
    //                         if (y->full) // A Hit for a matrix of length <= q
    //                         {
    //                                 if (hits[y->matrix] < max_hits){
    //                                 ret[y->matrix].push_back(match{i,y->score});
    //                                 hits[y->matrix]++;
    //                                 if (hits[y->matrix] == max_hits){
    //                                     matrices_left--;
    //                                     if (matrices_left == 0){
    //                                         return ret;
    //                                     }
    //                                 }
    //                             }

    //                             continue;
    //                         }
    //                         if (i - start >= motifs[y->matrix].window_pos() && i + motifs[y->matrix].size() - motifs[y->matrix].window_pos() <= end) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
    //                         {
    //                             double score;
    //                             bool m;

    //                             if (hits[y->matrix] < max_hits){

    //                                 std::tie(m,score) = motifs[y->matrix].check_hit(s, alphabet_map, i, y->score);

    //                                 if (m){
    //                                     ret[y->matrix].push_back(match{i - motifs[y->matrix].window_pos(),score});
    //                                     hits[y->matrix]++;
    //                                     if (hits[y->matrix] == max_hits){
    //                                         matrices_left--;
    //                                         if (matrices_left == 0){
    //                                             return ret;
    //                                         }
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
            
    //             // possible hits for matrices shorter than l near the end of current interval
    //             for (size_t i = end - l + 1; i < end; ++i)
    //             {
    //                 code = (code << SHIFT) & MASK;  // dummy character to the end of code
            
    //                 if (!window_hits[code].empty())
    //                 {
                
    //                     for (auto y = window_hits[code].begin(); y != window_hits[code].end(); ++y)
    //                     {
    //                         if (y->full && motifs[y->matrix].size() < end - i) // only sufficiently short hits are considered
    //                         {
    //                             if (hits[y->matrix] < max_hits){
    //                                 ret[y->matrix].push_back(match{i,y->score});
    //                                 hits[y->matrix]++;
    //                                 if (hits[y->matrix] == max_hits){
    //                                     matrices_left--;
    //                                     if (matrices_left == 0){
    //                                         return ret;
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         } 
    //     }
        
    //     return ret;
    // }

} // namespace scan
} // namespace MOODS
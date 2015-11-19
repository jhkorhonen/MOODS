// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include <utility>
#include <tuple>
#include <memory>
#include <algorithm>


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

        initialised = false;
        max_motif_size = 0;
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

        initialised = false;
        max_motif_size = 0;
    }
    
    void Scanner::set_motifs(const std::vector<score_matrix>& matrices,
                        const std::vector<double>& bg,
                        const std::vector<double> thresholds){

        this->motifs = vector<unique_ptr<Motif>>();
        
        for (size_t i = 0; i < matrices.size(); ++i){
            if (matrices[i].size() == a){
               motifs.emplace_back(new Motif0(matrices[i], bg, l, thresholds[i]));
            }
            else {
                motifs.emplace_back(new MotifH(matrices[i], bg, l, thresholds[i], a));
            }
            this->max_motif_size = std::max(this->max_motif_size, motifs.back()->size());
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

    
    // the template for various scan functions
    // the match_handler handles the actual logic of processing matches
    template<typename T> 
    void Scanner::process_matches(const std::string& s, T& match_handler)
    {

        if (!initialised){
            return;
        }

        const bits_t SHIFT = MOODS::misc::shift(a);
        const bits_t MASK = (1 << (SHIFT * l)) - 1;
        
        vector<size_t> bounds = misc::preprocess_seq(s, this->a, this->alphabet_map);
        
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
                    if (match_handler.has_hits(code))
                    {
                        code = (code << SHIFT) & MASK;  // dummy character to the end of code
                        
                        for (const scanner_output& y : match_handler.hits(code))
                        {
                            if (y.full && motifs[y.matrix]->size() <= end - i) // only sufficiently short hits are considered
                            {
                                match_handler.add_match(y.matrix, i, y.score);
                                // ret[y.matrix].push_back(match{i,y.score});
                            }
                        }
                        match_handler.clean_up();
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

                    if (match_handler.has_hits(code))
                    {
                        for (const scanner_output& y : match_handler.hits(code))
                        {
                            if (y.full) // A Hit for a matrix of length <= q
                            {
                                match_handler.add_match(y.matrix, i, y.score);
                                continue;
                            }
                            // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                            if (i - start >= motifs[y.matrix]->window_pos() && i + motifs[y.matrix]->size() - motifs[y.matrix]->window_pos() <= end)
                            {
                                double score;
                                bool m;
                                std::tie(m,score) = motifs[y.matrix]->check_hit(s, alphabet_map, i, y.score);
                                if (m){
                                    match_handler.add_match(y.matrix, i - motifs[y.matrix]->window_pos(), score);
                                }
                            }
                        }
                        match_handler.clean_up();
                    }
                }
            
                // possible hits for matrices shorter than l near the end of current interval
                for (size_t i = end - l + 1; i < end; ++i)
                {
                    code = (code << SHIFT) & MASK;  // dummy character to the end of code
            
                    if (match_handler.has_hits(code))
                    {
                        
                        for (const scanner_output& y : match_handler.hits(code))
                        {
                            if (y.full && motifs[y.matrix]->size() < end - i) // only sufficiently short hits are considered
                            {
                                match_handler.add_match(y.matrix, i, y.score);
                                // ret[y.matrix].push_back(match{i,y.score});
                            }
                        }
                        match_handler.clean_up();
                    }
                }
            } 
        }
    }


    // basic match handler
    // stores matches and returns all of them 
    class AllHitsMH{
    public:
        AllHitsMH(size_t motifs, std::vector<std::vector<scanner_output>>& _hits) : window_hits(_hits)
           {
            results = vector<vector<match> >(motifs, vector<match>());
           }

        bool has_hits(bits_t code){
            return !window_hits[code].empty();
        }

        vector<scanner_output>& hits(bits_t code){
            return window_hits[code];
        }

        void add_match(size_t matrix, size_t pos, double score){
            results[matrix].push_back(match{pos,score});
        }

        void clean_up(){}

        vector<vector<match>> get_results(){
            return results;
        }

    private:
        std::vector<std::vector<scanner_output>>& window_hits;
        std::vector<std::vector<match>> results;
    };

    std::vector<std::vector<match> > Scanner::scan(const std::string& s){
        AllHitsMH match_handler(motifs.size(), window_hits);
        process_matches<AllHitsMH> (s, match_handler);
        return match_handler.get_results();
    }

    class MaxHitsMH{
    public:
        MaxHitsMH(size_t motifs, std::vector<std::vector<scanner_output>>& _hits, size_t _max_hits)
           {
            window_hits = _hits; // copy?
            results = vector<vector<match> >(motifs, vector<match>());
            max_hits = _max_hits;

            to_clean = vector<size_t>();
            needs_cleaning = false;
           }

        bool has_hits(bits_t code){
            return !window_hits[code].empty();
        }

        vector<scanner_output>& hits(bits_t code){
            return window_hits[code];
        }

        void add_match(size_t matrix, size_t pos, double score){
            results[matrix].push_back(match{pos,score});
            if (results[matrix].size() >= max_hits){
                needs_cleaning = true;
                to_clean.push_back(matrix);
            }
        }

        void clean_up(){
            if (needs_cleaning){
                for (size_t code = 0; code < window_hits.size(); ++code){
                    for (size_t i = 0; i < to_clean.size(); ++i){
                        size_t matrix = to_clean[i];
                        for (auto y = window_hits[code].begin(); y < window_hits[code].end(); ++y){
                            if (y->matrix == matrix){
                                    window_hits[code].erase(y);
                                    break;
                            }
                        }
                    }
                }
                needs_cleaning = false;
                to_clean = vector<size_t>();
            }

        }

        vector<vector<match>> get_results(){
            return results;
        }

    private:
        std::vector<std::vector<scanner_output>> window_hits;
        std::vector<std::vector<match>> results;
        bool needs_cleaning;
        std::vector<size_t> to_clean;
        size_t max_hits;
    };

    std::vector<std::vector<match> > Scanner::scan_max_hits(const std::string& s, size_t max_hits){
        MaxHitsMH match_handler(motifs.size(), window_hits, max_hits);
        process_matches<MaxHitsMH> (s, match_handler);
        return match_handler.get_results();
    }

    class CountMaxHitsMH{
    public:
        CountMaxHitsMH(size_t motifs, std::vector<std::vector<scanner_output>>& _hits, size_t _max_hits)
           {
            window_hits = _hits; // copy?
            results = vector<size_t>(motifs, 0);
            max_hits = _max_hits;

            to_clean = vector<size_t>();
            needs_cleaning = false;
           }

        bool has_hits(bits_t code){
            return !window_hits[code].empty();
        }

        vector<scanner_output>& hits(bits_t code){
            return window_hits[code];
        }

        void add_match(size_t matrix, size_t pos, double score){
            results[matrix]++;
            if (results[matrix] >= max_hits){
                needs_cleaning = true;
                to_clean.push_back(matrix);
            }
        }

        void clean_up(){
            if (needs_cleaning){
                for (size_t code = 0; code < window_hits.size(); ++code){
                    for (size_t i = 0; i < to_clean.size(); ++i){
                        size_t matrix = to_clean[i];
                        for (auto y = window_hits[code].begin(); y < window_hits[code].end(); ++y){
                            if (y->matrix == matrix){
                                    window_hits[code].erase(y);
                                    break;
                            }
                        }
                    }
                }
                needs_cleaning = false;
                to_clean = vector<size_t>();
            }

        }

        vector<size_t> get_results(){
            return results;
        }

    private:
        std::vector<std::vector<scanner_output>> window_hits;
        std::vector<size_t> results;
        bool needs_cleaning;
        std::vector<size_t> to_clean;
        size_t max_hits;
    };

    std::vector<size_t> Scanner::counts_max_hits(const std::string& s, size_t max_hits){
        CountMaxHitsMH match_handler(motifs.size(), window_hits, max_hits);
        process_matches<CountMaxHitsMH> (s, match_handler);
        return match_handler.get_results();
    }


    void Scanner::variant_matches_recursive(std::vector<std::vector<match_with_variant>>& results, const state& current,
                                            const std::string& seq, const std::vector<variant> variants){


        size_t first_required_index = current.variant_start_pos + variants[current.vs.front()].modified_seq.size() - 1;         
        size_t last_required_index = current.prefix.size() - variants[current.vs.back()].modified_seq.size();

        size_t remaining = 0;
        if (current.prefix.size() - first_required_index < max_motif_size)
            remaining = max_motif_size - (current.prefix.size() - first_required_index);

        // first, we get modified hits for the current set of active sequence variants

        string current_seq = current.prefix + seq.substr(variants[current.vs.back()].end_pos, remaining);

        vector<vector<match>> current_results = this->scan(current_seq);
        
        for (size_t motif = 0; motif < current_results.size(); ++motif){
            for (size_t i = 0; i < current_results[motif].size(); ++i){
                match m = current_results[motif][i];

                if (m.pos + this->motifs[motif]->size() - 1 >= last_required_index && m.pos <= first_required_index) {
                    results[motif].push_back(match_with_variant{m.pos + current.seq_start_pos,m.score,current.vs});
                }
            }    
        }

        // second, recursive call for each possible variant that can follow the current last variant

        for (size_t next_variant = current.vs.back() + 1; next_variant < variants.size(); ++next_variant){

            variant cv = variants[current.vs.back()];
            variant nv = variants[next_variant];

            if (nv.start_pos - cv.end_pos >= remaining){
                break;
            }

            // add next variant and recurse if the next variant is compatible with the current one AND
            // it is not a deletion at one end of the sequence
            if (
                 ((cv.end_pos < nv.start_pos) ||
                  (cv.end_pos == nv.start_pos && ( cv.start_pos < cv.end_pos || nv.start_pos < nv.end_pos  )))
                   &&
                !((nv.start_pos == 0 && nv.modified_seq.size() == 0) ||
                  (nv.end_pos == seq.size() && nv.modified_seq.size() == 0))
                )
            {

                string next_prefix = current.prefix + seq.substr(cv.end_pos, nv.start_pos - cv.end_pos) + nv.modified_seq;
                auto next_variants = current.vs;
                next_variants.push_back(next_variant);

                state next = {next_variants, next_prefix, current.seq_start_pos, current.variant_start_pos};

                this->variant_matches_recursive(results, next, seq, variants);
            }
        }        
    }


    std::vector<std::vector<match_with_variant>> Scanner::variant_matches(const std::string& seq, const std::vector<variant> _variants){
        
        // copy and sort the variants
        std::vector<variant> variants = _variants;
        std::sort(variants.begin(), variants.end());

        vector<vector<match_with_variant>> results(this->size(), vector<match_with_variant>());

        for (size_t i = 0; i < variants.size(); ++i){
            // we'll skip all deletions that start at position 0 or 
            // end at the last position 
            // as that does not really make any sense in our framework
            if (! ( (variants[i].start_pos == 0 && variants[i].modified_seq.size() == 0) ||
                    (variants[i].end_pos == seq.size() && variants[i].modified_seq.size() == 0))) {
                        size_t prefix_start = 0;
                        size_t prefix_length = 0;
                        if (variants[i].start_pos >= max_motif_size-1){
                            prefix_start = variants[i].start_pos - max_motif_size + 1;
                            prefix_length = max_motif_size - 1;
                        }
                        else {
                            prefix_start = 0;
                            prefix_length = variants[i].start_pos;
                        }
        
                        string next_prefix = seq.substr(prefix_start, prefix_length) + variants[i].modified_seq;
        
                        state s = {vector<size_t>(1,i), // active variants
                                   next_prefix, // prefix string 
                                   prefix_start, // seq position where prefix starts
                                   prefix_length // first position in prefix that is from the modified string 
                                  };
        
                        this->variant_matches_recursive(results, s, seq, variants);
            }
        }

        return results;
    }
    



} // namespace scan
} // namespace MOODS
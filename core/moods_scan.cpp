

#include "moods.h"
#include "moods_scan.h"
#include "moods_misc.h"
#include "moods_tools.h"
#include "match_types.h"
#include "motif.h"
#include "scanner.h"

using std::vector;
using std::string;
using std::size_t;


namespace MOODS { namespace scan{

    vector< vector<match> > scan_dna(const string& seq, const vector<score_matrix>& matrices,
                                            const vector<double>& bg, const vector<double> thresholds, unsigned int window_size )
    {
        
        Scanner scanner(window_size);
        scanner.set_motifs(matrices, bg, thresholds);

        auto results = scanner.scan(seq);        
        return results;
    }

    std::vector< std::vector< match> > scan(const std::string& seq, const std::vector<score_matrix>& matrices,
                                            const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size,
                                            const std::vector<std::string>& alphabet)
    {

        Scanner scanner(window_size, alphabet);
        scanner.set_motifs(matrices, bg, thresholds);
        
        auto results = scanner.scan(seq);
        return results;
    }

    std::vector< std::vector< match> > scan_best_hits_dna(const std::string& seq, const std::vector<score_matrix>& matrices, size_t target,
                                                          int iterations, unsigned int MULT, size_t LIMIT_MULT, size_t window_size){

        vector<double> bg = tools::bg_from_sequence_dna(seq,0.01);

        // try to guess a good initial  threshold; this works quite often
        double p = ((1 + 2.0*MULT) / 3.0) * target/ seq.size();

        if (LIMIT_MULT < MULT){
            MULT = LIMIT_MULT;
        }

        vector<double> current (matrices.size(),0);
        vector<double> upper (matrices.size(),0);
        vector<double> lower (matrices.size(),0);

        vector<bool> ok(matrices.size(),0);
        size_t remaining = matrices.size();

        for (size_t i = 0; i < matrices.size(); ++i){
            current[i] = tools::threshold_from_p(matrices[i],bg,p);
            upper[i] = tools::max_score(matrices[i],4);
            lower[i] = tools::min_score(matrices[i],4);
        }

        bool done = 0;

        int iteration = 0;


        // binary search for all matrices in parallel
        while (remaining > 0){
            iteration++;

            vector<size_t> it_indices;
            vector<score_matrix> it_matrices;
            vector<double> it_thresholds;
            vector<char> status (matrices.size(), '_');

            for (size_t i = 0; i < matrices.size(); ++i){
                if (!ok[i]){
                    it_indices.push_back(i);
                    it_matrices.push_back(matrices[i]);
                    it_thresholds.push_back(current[i]);
                }
            }

            Scanner scanner(window_size);
            scanner.set_motifs(it_matrices, bg, it_thresholds);
        
            vector<size_t> results = scanner.counts_max_hits(seq, LIMIT_MULT * target);

            for (size_t j = 0; j < it_indices.size(); ++j){
                size_t i = it_indices[j];
                size_t hits = results[j];
                double threshold = current[i];
                if (hits >= target && hits < MULT*target){
                    ok[i] = 1;
                    remaining--;
                }
                // too many hits, INCREASE the threshold & lower bound
                else if (hits >= MULT*target){
                    status[i] = 'H';
                    current[i] = (upper[i] + threshold)/2.0;
                    lower[i] = threshold;
                }
                // too few this, DECREASE the threshold & upper bound
                else{
                    status[i] = 'L';
                    current[i] = (lower[i] + threshold)/2.0;
                    upper[i] = threshold;
                }
                // failsafe
                if ((current[i] == threshold || iteration == iterations) && !ok[i]){
                    ok[i] = 1;
                    status[i] = 'X';
                    if (hits >= LIMIT_MULT * target){
                        current[i] = upper[i];
                    }
                    remaining--;
                }
            }
        }



        Scanner scanner(window_size);
        scanner.set_motifs(matrices, bg, current);
        
        return scanner.scan(seq); 

    }

    std::vector<match> naive_scan_dna(const std::string& seq, const score_matrix matrix, double threshold){

        vector<unsigned char> alphabet_map (256, 4);
        
        alphabet_map[(unsigned char)'a'] = 0;
        alphabet_map[(unsigned char)'A'] = 0;
        
        alphabet_map[(unsigned char)'c'] = 1;
        alphabet_map[(unsigned char)'C'] = 1;
        
        alphabet_map[(unsigned char)'g'] = 2;
        alphabet_map[(unsigned char)'G'] = 2;
        
        alphabet_map[(unsigned char)'t'] = 3;
        alphabet_map[(unsigned char)'T'] = 3; 

        size_t m = matrix[0].size();
        vector<match> results;

        vector<size_t> bounds = misc::preprocess_seq(seq, 4, alphabet_map);
        
        // Scanning
        for (size_t seq_i = 0; seq_i < bounds.size(); ){
            size_t start = bounds[seq_i];
            ++seq_i;
            size_t end = bounds[seq_i];
            ++seq_i;

            for (size_t i = start; i + m < end + 1; ++i){
                double score = 0;

                for (size_t j = 0; j < m; ++j){
                    score += matrix[alphabet_map[seq[i+j]]][j];
                }

                if (score >= threshold){
                    results.push_back(match{i,score});
                }

            }
        }

        return results;
    }

    std::vector<match> naive_scan_dna(const std::string& seq, const score_matrix matrix, double threshold, size_t a){

        vector<unsigned char> alphabet_map (256, 4);
        
        alphabet_map[(unsigned char)'a'] = 0;
        alphabet_map[(unsigned char)'A'] = 0;
        
        alphabet_map[(unsigned char)'c'] = 1;
        alphabet_map[(unsigned char)'C'] = 1;
        
        alphabet_map[(unsigned char)'g'] = 2;
        alphabet_map[(unsigned char)'G'] = 2;
        
        alphabet_map[(unsigned char)'t'] = 3;
        alphabet_map[(unsigned char)'T'] = 3; 

        size_t cols = matrix[0].size();
        size_t rows = matrix.size();
        vector<match> results;

        size_t q = MOODS::misc::q_gram_size(rows, a);
        bits_t SHIFT = MOODS::misc::shift(a);
        bits_t MASK = (1 << (SHIFT * q)) - 1;

        vector<size_t> bounds = misc::preprocess_seq(seq, 4, alphabet_map);
        
        // Scanning
        for (size_t seq_i = 0; seq_i < bounds.size(); ){
            size_t start = bounds[seq_i];
            ++seq_i;
            size_t end = bounds[seq_i];
            ++seq_i;

            for (size_t i = start; i + cols + q - 1 < end; ++i){
                double score = 0;
                bits_t CODE = 0;

                for (size_t k = 0; k < q - 1; ++k){
                    CODE = (CODE << SHIFT) | alphabet_map[seq[i+k]];
                }

                for (size_t j = 0; j < cols; ++j){
                    CODE = MASK & (CODE << SHIFT) | alphabet_map[seq[i+j+q-1]];
                    score += matrix[CODE][j];
                }

                if (score >= threshold){
                    results.push_back(match{i,score});
                }

            }
        }

        return results;
    }
}}


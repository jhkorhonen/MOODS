

#include "moods.h"
#include "moods_scan.h"
#include "moods_misc.h"
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


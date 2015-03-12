
#include <iostream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <cstdlib>
#include <math.h>


#include "pssm_algorithms.hpp"

using std::cerr;
using std::cout;

using std::cerr;
using std::cout;

// Struct for comparing matrix rows with sort()
struct compareRows
{
    const doubleArray *goodness;
    bool operator() (int i, int j)
    {
        return ( (*goodness)[i] > (*goodness)[j] );
    }
};

struct OutputListElementMulti
{
    score_t score;
    int matrix;
    bool full;
};


// Transforms a weight matrix into a PSSM
scoreMatrix freq2LogOddsHO(const scoreMatrix &mat, const doubleArray &bg)
{
    int numA = 16;
    int n = mat[0].size();

    std::vector<score_t> col(n, 0);
    scoreMatrix ret(numA, col);
    
    
    // 1-order terms
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < numA; ++j)
        {
            ret[j][i] = log(mat[j][i]) - log(bg[j & 3]);
        }
    }
    
    // 0-order terms to the first row
    for (int j = 0; j < 4; ++j)
    {
        for (int k = 0; k < 4; ++k)
        {
            ret[(j << 2) & k][0] += log(mat[numA & j][0]) - log(bg[j]);
        }
    }
    
    return ret;
}

scoreArray prefix_scores(const scoreMatrix &mat, const int q, const int end){
    // the "actual" strings are 0...(end-1)
    // the corresponding indices of the matrix are 0...(end-q)
    //
    // 
    // 012345678
    // |--| <- end = 4, q = 2
    // 00224466
    //  11335577
    
    const int numA = 1 << (2*q); // 4**q
    const int numQ = 1 << (2*(q-1));
    const int maskQ = numQ - 1;
    
    scoreArray p_scores(numQ, 0);
    
    // we should have end >= q
    // otherwise we assume that we want 0 as an aswer
    if (end < q){
        return p_scores; // all zeroes
    }
    
    for (int j = 0; j < end-q; ++j){
        scoreArray p_scores_n(numQ, SCORE_MIN);
        
        for (int i = 0; i < numA; ++i){
            p_scores_n[i & maskQ] = std::max(mat[i][j] + p_scores[i >> 2], p_scores_n[i & maskQ]);
        }
        
        p_scores = p_scores_n;
    }
    
    return p_scores;
}

scoreArray suffix_scores(const scoreMatrix &mat, const int q, const int start){
    // the "actual" strings are start...(m+q-1)
    // the corresponding indices of the matrix are start...(m-1)
    //
    // 
    // 012345678
    //      |--| <- start = 5, q = 2
    // 00224466
    //  11335577
    
    const int numA = 1 << (2*q); // 4**q
    const int numQ = 1 << (2*(q-1));
    const int maskQ = numQ - 1;
    
    const int m = mat[0].size();
    
    scoreArray s_scores(numQ, 0);
    
    if (start > m-1){
        return s_scores; // all zeroes
    }
    
    for (int j = m-1; j >= start; --j){
        scoreArray s_scores_n(numQ, SCORE_MIN);
        
        for (int i = 0; i < numA; ++i){
            s_scores_n[i >> 2] = std::max(mat[i][j] + s_scores[i & maskQ], s_scores_n[i >> 2]);
        }
        
        s_scores = s_scores_n;
    }
    
    return s_scores;
}

scoreArray window_scores(const scoreMatrix &mat, const int q, const int l){
    const int m = mat[0].size(); // actual length - q + 1
    const int numA = 1 << (2*q); // 4**q
    const int numQ = 1 << (2*(q-1));
    const int maskQ = numQ - 1;
    
    scoreArray ret(m - l + 1, SCORE_MIN);
    
    for (int start = 0; start < m - l-1; ++start){
        scoreArray p_scores(numQ, SCORE_MIN);
        for (int i = 0; i < numA; ++i){
            p_scores[i & maskQ] = std::max(mat[i][start], p_scores[i & maskQ]);
        }
        
        for (int j = start + 1; j < start + l - q - 2; ++j){
            scoreArray p_scores_n(numQ, SCORE_MIN);
        
            for (int i = 0; i < numA; ++i){
                p_scores_n[i & maskQ] = std::max(mat[i][j] + p_scores[i >> 2], p_scores_n[i & maskQ]);
            }
            p_scores = p_scores_n;
        }
        
        for (int i = 0; i < numQ; ++i){
            ret[start] = std::max(p_scores[i], ret[start]);
        }
    }
    
    return ret;

}

scoreMatrix expected_scores(const std::vector<scoreMatrix> &matrices, const doubleArray &bg, const int numA){
    scoreMatrix goodnesses;
    goodnesses.reserve(matrices.size());
    
    intArray m(matrices.size(), 0);
    for (int i = 0; i < (int) matrices.size(); ++i)
    {
        m[i] = matrices[i][0].size();
    }
    
    
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        doubleArray expd(m[k]);

        for (int i = 0; i < m[k]; ++i)
        {
            for (int j = 0; j < numA; ++j)
            {
                expd[i] = bg[j ^ 3] * bg[(j >> 2) & 3]  * matrices[k][j][i];
            }
        }
        goodnesses.push_back(expd);
    }
    
    return goodnesses;
}

scoreMatrix expected_losses(const std::vector<scoreMatrix> &matrices, const doubleArray &bg, const int numA){

    scoreMatrix exp_loss = expected_scores(matrices, bg, numA); // get the expected scores
    
    for (int k = 0; k < (int)matrices.size(); ++k){
        for (int i = 0; i < (int)exp_loss[k].size(); ++i){
            score_t max = SCORE_MIN;
            for (int j = 0; j < numA; ++j)
            {
                max = std::max(max, matrices[k][j][i]);
            }
            
            exp_loss[k][i] = max - exp_loss[k][i];
        }
    }
    
    return exp_loss;
}

// Multimatrix LFA for DNA
// vector<matchArray> mmlf_DNA_HO(const int q, const charArray &s, const vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol)
std::vector<matchArray> mmlf_DNA_HO(const int q, const charArray &s, const std::vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol)
{

    const int BITSHIFT = 2;
    const unsigned int numA = 4; // 2**BITSIFT

    const position_t n = s.size();

    intArray m(matrices.size(), 0);
    for (int i = 0; i < (int) matrices.size(); ++i)
    {
        m[i] = matrices[i][0].size();
    }

    // Calculate expected differences
    scoreMatrix exp_score = expected_scores(matrices,bg,numA);
    
    intArray window_positions;
    window_positions.reserve(matrices.size());
    
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            window_positions.push_back(0);
        }
        else
        {
            scoreArray w_scores = window_scores(matrices[k],2,q);
            double current_exp = 0;
            for (int i = 0; i < q-1; ++i)
            {
                current_exp += exp_score[k][i];
            }

            double max_loss = w_scores[0] - current_exp;
            int window_pos = 0;
            
            cerr << w_scores[0] - current_exp << "\n";

            for (int i = 0; i < m[k] - q + 1; ++i)
            {
                current_exp -= exp_score[k][i];
                current_exp += exp_score[k][i+q-1];
                
                cerr << w_scores[i+1] - current_exp << "\n";
                if (w_scores[i+1] - current_exp > max_loss)
                {
                    max_loss = w_scores[i+1] - current_exp;
                    window_pos = i+1;
                }
            }
            window_positions.push_back(window_pos);
            cerr << k << " " << window_pos << "\n";
        }
    }
    
    // Pre- and post-window scores
    scoreMatrix P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        scoreArray B = prefix_scores(matrices[k],2,window_positions[k]+2);
        P.push_back(B);
    }
    
    scoreMatrix S;
    S.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        scoreArray B = suffix_scores(matrices[k],2,window_positions[k]+q-2);
        S.push_back(B);
    }
    
    // Arrange matrix indices not in window by expected loss, for use in scanning
    intMatrix orders;
    orders.reserve(matrices.size());
    scoreMatrix L;
    L.reserve(matrices.size());
    
    scoreMatrix exp_losses = expected_losses(matrices,bg,numA);

    for (unsigned short k = 0; k < (int) matrices.size(); ++k)
    {
        if (q >= m[k]+1)
        {
            intArray temp1;
            orders.push_back(temp1);
            scoreArray temp2;
            L.push_back(temp2);
        }
        else
        {
            for (int i = 0; i < exp_losses[k].size(); ++i){
                cerr << i << " " << exp_losses[k][i] << "\n";
            }
            intArray order(m[k]-q+1, 0);
            for (int i = 0; i < window_positions[k]; ++i)
            {   
                order[i] = i;
            }
            for (int i = window_positions[k]+q-1; i < m[k]; ++i)
            {
                order[i-q+1] = i;
            }

            compareRows comp;
            comp.goodness = &(exp_losses[k]);

            sort(order.begin(), order.end(), comp);

            orders.push_back(order);
            
            for (int i = 0; i < order.size();++i){
                cerr << order[i] << " ";
            }
            cerr << "\n";

            scoreArray K(m[k]-q+1, 0);
            for (int j = m[k]-q; j > 0; --j)
            {
                score_t max = INT_MIN;
                for (unsigned int i = 0; i < numA; ++i)
                {
                    max = std::max(max,matrices[k][i][order[j]]);
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }
    
    //
    const bits_t size = 1 << (BITSHIFT * q); // 4^q -- the actual size of the window
    const bits_t BITAND = size - 1;
    
    std::vector<std::vector< OutputListElementMulti> > output(size);

    {
        bitArray sA(q,0);
        while (true)
        {
            bits_t code = 0;
            for (int j = 0; j < q; ++j)
            {
                code = (code << BITSHIFT) | sA[j];
            }

            for (unsigned int k = 0; k < matrices.size(); ++k )
            {
                if (m[k] + 1 <= q)
                {
                    score_t score = 0;
                    for (int i = 0; i < m[k]; ++i)
                    {
                        score += matrices[k][(sA[i] << 2) | sA[i+1]][i];
                    }
                    if (score >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = true;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
                else
                {
                    score_t score = 0;
                    for (int i = 0; i < q-1; ++i)
                    {
                        score += matrices[k][(sA[i] << 2) | sA[i+1]][i + window_positions[k]];
                    }
                    if (score + P[k][sA[0]] + S[k][sA[q-1]] >= tol[k])
                    {
                        OutputListElementMulti temp;
                        temp.full = false;
                        temp.matrix = k;
                        temp.score = score;
                        output[code].push_back(temp);
                    }
                }
            }

            int pos = 0;
            while (pos < q)
            {
                if (sA[pos] < numA - 1)
                {
                    ++sA[pos];
                    break;
                }
                else
                {
                    sA[pos] = 0;
                    ++pos;
                }
            }

            if (pos == q)
                break;
        }
    }
    //
    // // Scanning
    //
        
    std::vector<matchArray> ret;
        for (unsigned int i = 0; i < matrices.size(); ++i)
        {
            matchArray temp;
            ret.push_back(temp);
        }

    matchData hit;

    score_t score;
    position_t k;
    position_t limit;
    position_t ii;
    score_t tolerance;
    intArray::iterator z;

    bits_t code = 0;
    for (position_t ii = 0; ii < q - 1; ++ii)
        code = (code << BITSHIFT) + s[ii];


    for (position_t i = 0; i < n - q + 1; ++i)
    {
        code = ((code << BITSHIFT) + s[i + q - 1]) & BITAND;

        if (!output[code].empty())
        {
            for (auto y = output[code].begin(); y != output[code].end(); ++y)
            {
                if (y->full) // A Hit for a matrix of length <= q
                {
                    hit.position = i;
                    hit.score = y->score;
                    ret[y->matrix].push_back(hit);
                    continue;
                }
                if (i - window_positions[y->matrix] >= 0 && i + m[y->matrix] + 1 - window_positions[y->matrix] <= n) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                {
                    score = y->score;
                    k = y->matrix;
                    limit = m[k] - q + 1;
                    ii = i - window_positions[k]; 
                    tolerance = tol[k];
                    z = orders[k].begin();
                    for (int j = 0; j < limit  ;++j)
                    {
                        score += matrices[k][(s[ii+(*z)] << 2) | s[ii+(*z) + 1]][*z];
                        // if (score + L[k][j] < tolerance)
                        //     break;
                        ++z;
                    }
                    if (score >= tolerance){
                    hit.position = i - window_positions[k];
                    hit.score = score;
                    ret[k].push_back(hit);
                }
                }
            }
        }
    }

    for (position_t i = n - q + 1; i < n; ++i) // possible hits for matrices shorter than q near the end of sequence
    {
        code = (code << BITSHIFT) & BITAND; // dummy character to the end of code

        if (!output[code].empty())
            {

            for (auto y = output[code].begin(); y != output[code].end(); ++y)
                {
                    if (y->full && m[y->matrix] + 1 < n - i + 1) // only sufficiently short hits are considered
                    {
                        hit.position = i;
                        hit.score = y->score;
                        ret[y->matrix].push_back(hit);
                    }
                }
        }
    }

    return ret;
}

matchArray naive_HO(const charArray &s, const scoreMatrix &p, const score_t tol)
{
    const int m = p[0].size();
    const position_t n = s.size();

    matchArray ret;
    matchData hit;

    for (position_t i = 0; i <= n - m; ++i)
    {
        score_t score = 0;
        for (int j = 0; j < m; ++j)
        {
            score += p[ (s[i + j] << 2) | s[i + j + 1] ][j];
        }
        if (score >= tol)
        {
            hit.position = i;
            hit.score = score;
            ret.push_back(hit);
        }
    }

    return ret;
}


int main(int argc, char **argv)
{
    int c;
    int q = 4;
    double tol;
    double pseudocount = 1;
    
    bool print_results = true;
    clock_t single;
    clock_t total;

    doubleArray bg;

    while (1)
        {
        static struct option long_options[] =
            {
                {"help",       no_argument,       0, 'h'},
                {"parameter",  required_argument, 0, 'p'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "a:hp:fqb:g:",
                            long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
            {
            case 'h':
                cerr << "usage: pml [options] p-value sequence matrix1 matrix2 ...\n";
                exit(0);
                break;
            case 'p':
                q = atoi(optarg);
                if (q < 1)
                {
                    cerr << "Invalid parameter p: " << q << "\n";
                    exit(1);
                }
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                exit(1);
            }
        }


    if (!(optind + 2 < argc))
    {
        cerr << "usage: test_hom [options] p-value sequence matrix1 matrix2 ...\n";
        exit(1);
    }

    tol = atof(argv[optind]);
    
    ++optind;

    std::ifstream ifs(argv[optind]);
    if(!ifs) {
    	cerr <<"Could not open file "<<argv[optind]<<"\n";
    	exit(1);
    }
    
    auto seq = readDNA(ifs);

    ++optind;



    bg = bgFromSequence(seq, 4, pseudocount);

    cerr << "calculated background\n";

    std::vector<scoreMatrix> matrices;
    std::vector<score_t> thresholds;
    std::vector< char** > names;

    for (int i = optind; i < argc; ++i)
    {
        std::ifstream ifs(argv[i]);
        if (!ifs)
        {
            cerr << "Error reading matrices\n";
            exit(1);
        }
        scoreMatrix matrix = readMatrix(ifs);
        matrices.push_back(  freq2LogOddsHO( matrix, bg)  );
        thresholds.push_back(tol);
        names.push_back(&(argv[i]));
        cerr << "loaded matrix " << argv[i] << "\n";
        printScoreMatrix(matrix);
        printScoreMatrix(matrices.back());
    }

    cerr << "loaded matrices\n";

    total = clock();
    position_t total_results = 0;

    cerr << "Naive algorithm\n";
    for (unsigned int k = 0; k < matrices.size(); ++k){
        cout << *(names[k]) << "\n";
        matchArray results;
        results = naive_HO(seq, matrices[k], thresholds[k]);
        cerr << "Hits: " << results.size() << "\n";
        printMatchArray(results);
    }
    
    
    std::vector<matchArray> results;
    results = mmlf_DNA_HO(q, seq, matrices, bg, thresholds);
    
    for (unsigned int k = 0; k < matrices.size(); ++k)
    {
        cerr << *(names[k]) << "\n";
        cerr << "Hits: " << results[k].size() << "\n";
        printMatchArray(results[k]);
        cout << "\n";
    }


        //     if(algorithm < 4)
        //     {
        //         for (unsigned int i = 0; i < matrices.size(); ++i)
        //         {
        //             cout << *(names[i]) << "\n";
        //             matchArray results;
        //
        //             seq_it.reset();
        //             single = clock();
        //             switch (algorithm)
        //             {
        //                 case 0: results = naiveAlgorithm(seq_it, matrices[i], thresholds[i]); break;
        //                 case 1: results = naiveSuperalphabetAlgorithmDNA(q, seq_it, matrices[i], thresholds[i]); break;
        //                 case 2: results = permutatedLookAhead(seq_it, matrices[i], bg, thresholds[i]); break;
        //                 case 3: results = lookaheadFiltrationDNA(q, seq_it, matrices[i], bg, thresholds[i]); break;
        //
        //             }
        //
        //             single = clock() - single;
        //
        //             cerr << "Time: " << std::setprecision(5) << (double)single/((double)CLOCKS_PER_SEC) << "\n";
        //             cerr << "Hits: " << results.size() << "\n";
        //             total_results += results.size();
        //             if(print_results)
        //                 printMatchArray(results);
        //             cout << "\n";
        //         }
        //     }
        //     else
        //     {
        //         std::vector<matchArray> results;
        //         switch (algorithm)
        //         {
        //             case 4: results = multipleMatrixLookaheadFiltrationDNA(q, seq_it, matrices, bg, thresholds); break;
        //             case 5: results = searchDNA(seq_it, matrices, thresholds, bg, q); break;
        //         }
        //
        // for (unsigned int i = 0; i < matrices.size(); ++i)
        // {
        //     cerr << *(names[i]) << "\n";
        //     cerr << "Hits: " << results[i].size() << "\n";
        //     total_results += results[i].size();
        //     if(print_results)
        //         printMatchArray(results[i]);
        //     cout << "\n";
        // }
        //     }

    total = clock() - total;

    // cerr << "Total time: " << std::setprecision(5) << (double)total/((double)CLOCKS_PER_SEC) << "\n";
    // cerr << "Total results: " << total_results <<"\n";
    // cout << "\n";


    return 0;
}


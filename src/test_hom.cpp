
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

// Multimatrix LFA for DNA
vector<matchArray> mmlf_DNA_HO(const int q, const charArray &s, const vector<scoreMatrix> &matrices, const doubleArray &bg, const scoreArray &tol)
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
    vector<doubleArray> goodnesses;
    goodnesses.reserve(matrices.size());

    for (int i = 0; i < (int)matrices.size(); ++i)
    {
        doubleArray expd(m[i]);

        for (int i = 0; i < m; ++i)
        {
            score_t max = SCORE_MIN;
            for (int j = 0; j < numA; ++j)
            {
                if (max < mat[j][i])
                    max = mat[j][i];
            }

            expd[i] = max;

            for (int j = 0; j < numA; ++j)
            {
                expd[i] -= bg[j ^ 3] * bg[(j >> 2) & 3]  * mat[j][i];
            }
        }
        goodnesses.push_back(expd);
    }
    
    // --- fix everything below this line ---
    
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
            double current_goodness = 0;
            for (int i = 0; i < q; ++i)
            {
                current_goodness += goodnesses[k][i];
            }

            double max_goodness = current_goodness;
            int window_pos = 0;

            for (int i = 0; i < m[k] - q; ++i)
            {
                current_goodness -= goodnesses[k][i];
                current_goodness += goodnesses[k][i+q];
                if (current_goodness > max_goodness)
                {
                    max_goodness = current_goodness;
                    window_pos = i+1;
                }
            }
            window_positions.push_back(window_pos);
        }
    }

    // Calculate lookahead scores for all matrices
    scoreMatrix T;
    T.reserve(matrices.size());

    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        scoreArray C(m[k],0);
        for (int j = m[k] - 1; j > 0; --j) {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T.push_back(C);
    }

    // Pre-window scores
    scoreArray P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        score_t B = 0;
        for (int j = 0; j < window_positions[k]; ++j)
        {
            score_t max = SCORE_MIN;
            for (unsigned int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            B += max;
        }
        P.push_back(B);
    }

    // Arrange matrix indeces not in window by entropy, for use in scanning
    intMatrix orders;
    orders.reserve(matrices.size());
    scoreMatrix L;
    L.reserve(matrices.size());

    for (unsigned short k = 0; k < (int) matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            intArray temp1;
            orders.push_back(temp1);
            scoreArray temp2;
            L.push_back(temp2);
        }
        else
        {
            intArray order(m[k]-q, 0);
            for (int i = 0; i < window_positions[k]; ++i)
            {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i)
            {
                order[i-q] = i;
            }

            compareRows comp;
            comp.goodness = &(goodnesses[k]);

            sort(order.begin(), order.end(), comp);

            orders.push_back(order);

            scoreArray K(m[k]-q, 0);
            for (int j = m[k]-q-1; j > 0; --j)
            {
                score_t max = INT_MIN;
                for (unsigned int i = 0; i < numA; ++i)
                {
                    if (max < matrices[k][i][order[j]])
                        max = matrices[k][i][order[j]];
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }

    const bits_t size = 1 << (BITSHIFT * q); // numA^q
    const bits_t BITAND = size - 1;


    vector<vector< OutputListElementMulti> > output(size);

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
                if (m[k] <= q)
                {
                    score_t score = 0;
                    for (int i = 0; i < m[k]; ++i)
                    {
                        score += matrices[k][sA[i]][i];
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
                    for (int i = 0; i < q; ++i)
                    {
                        score += matrices[k][sA[i]][i + window_positions[k]];
                    }
                    if (score + P[k] + T[k][q + window_positions[k]-1] >= tol[k])
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

	// Scanning

	vector<matchArray> ret;
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
            for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
            {
                if (y->full) // A Hit for a matrix of length <= q
                {
                    hit.position = i;
                    hit.score = y->score;
                    ret[y->matrix].push_back(hit);
                    continue;
                }
                if (i - window_positions[y->matrix] >= 0 && i + m[y->matrix] - window_positions[y->matrix] <= n) // A possible hit for a longer matrix. Don't check if matrix can't be positioned entirely on the sequence here
                {
                    score = y->score;
                    k = y->matrix;
                    limit = m[k] - q;
                    ii = i - window_positions[k];
                    tolerance = tol[k];
                    z = orders[k].begin();
                    for (int j = 0; j < limit  ;++j)
                    {
                        score += matrices[k][s[ii+(*z)]][*z];
                        if (score + L[k][j] < tolerance)
                            break;
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

			for (vector<OutputListElementMulti>::iterator y = output[code].begin(); y != output[code].end(); ++y)
        	{
            	if (y->full && m[y->matrix] < n - i + 1) // only sufficiently short hits are considered
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


int main(int argc, char **argv)
{
    int c;
    int q = 7;
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


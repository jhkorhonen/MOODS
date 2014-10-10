// find_pssm_dna - tool for finding PSSM matches from dna sequences
// Copyright (C) 2007-2009  Pasi Rastas, Janne Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.



#include <iostream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <cstdlib>

#include "pssm_algorithms.hpp"
#include "seq_buffer.h"

class SeqStreamSource : public SeqSourceI {
	char c;
	int value;
	std::istream &in;
	static const int FILE_BUFFER_SIZE = 1024;
	char buffer[FILE_BUFFER_SIZE];
public:
	SeqStreamSource(std::istream &input_stream) : in(input_stream) {
	}
	int read_data(char *dst, int length) {
		int i = 0;
		bool reading_head = false;
		while(i < length && (!in.eof())) {
			in.read(buffer, (length - i < FILE_BUFFER_SIZE) ? length - i : FILE_BUFFER_SIZE);
			if(in.gcount() == 0) {
				break;
			}
			int k = 0;
			//printf("stream read %d input %d\n", (length - i < FILE_BUFFER_SIZE) ? length - i : FILE_BUFFER_SIZE, in.gcount());
			for( ; k < in.gcount(); k++) {
				if(!reading_head) {
					switch (buffer[k])
					{
						case 'a':
						case 'A': dst[i] = 0; ++i;break;
						case 'c':
						case 'C': dst[i] = 1; ++i;break;
						case 'g':
						case 'G': dst[i] = 2; ++i; break;
						case 't':
						case 'T': dst[i] = 3; ++i; break;
						case 'n':
						case 'N': dst[i] = (int) (4 * rand() / (1.0 + RAND_MAX)); ++i; break;
						case '>':
							reading_head = true;
							break;//you should break for loop too
						default:
							break;
					}
				}
				else {
					if(buffer[k] == '\n')
					{
							  reading_head = false;
					}
				}
			}
		}
	    //printf("read_data called %d current pos %ld data read %d\n", length, (unsigned long) in.tellg(),  i);
	    return i;
	}
	void reset() {
		in.clear();
		in.seekg(std::ios::beg);
	}
	bool eof() {
		return in.eof();
	}
};

using std::cerr;
using std::cout;

int main(int argc, char **argv)
{
    int c;
    int algorithm = 5;
    int q = 7;
    double tol;
    bool flat_bg = false;
    bool print_results = true;
    int buffer_size = 60000000;
    clock_t single;
    clock_t total;

    doubleArray bg;

    while (1)
        {
        static struct option long_options[] =
            {
                {"help",       no_argument,       0, 'h'},
                {"flat_bg",    no_argument,       0, 'f'},
                {"quiet",      no_argument,       0, 'q'},
                {"algorithm",  required_argument, 0, 'a'},
                {"parameter",  required_argument, 0, 'p'},
                {"buffer",     required_argument, 0, 'b'},
                {"bg_file",	   required_argument, 0, 'g'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "a:hp:fqb:g:",
                            long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        std::ifstream bg_in;
        switch (c)
            {
            case 'h':
                cerr << "usage: pml [options] p-value sequence matrix1 matrix2 ...\n";
                exit(0);
                break;

            case 'a':
                {
                    int a = atoi(optarg);
                    if (a >= 0 && a < 5)
                        algorithm = a;
                    else
                    {
                        cerr << "Algorithms: 0=naive, 1=ns, 2=pls, 3=lf, 4=mlf\n";
                        exit(1);
                    }
                }
                break;

            case 'p':
                q = atoi(optarg);
				if (q < 1)
				{
					cerr << "Invalid parameter p: " << q << "\n";
					exit(1);
				}
                break;
            case 'b':
            	buffer_size = atoi(optarg);
				if (buffer_size < 1)
				{
					cerr << "Invalid parameter b: " << buffer_size << "\n";
					exit(1);
				}
                break;
            case 'f':
                flat_bg = true;
                break;
            case 'q':
            	print_results = false;
            	break;
            case 'g':
            	bg_in.open(optarg);
            	if(!bg_in) {
            		cerr<<"Could not open bg file "<<optarg<<std::endl;
            		exit(1);
            	}
            	for(int i = 0; i < 4; i++) {
					double tmp;
					bg_in>>tmp;
					if(tmp < 0 || tmp > 1) {
						cerr<<"Invalid background distribution file"<<std::endl;
						exit(1);
					}
					bg.push_back(tmp);
            	}
            	bg_in.close();
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
        cerr << "usage: pml [options] p-value sequence matrix1 matrix2 ...\n";
        exit(1);
    }

    tol = atof(argv[optind]);
	if(tol < 0 || tol > 1)
	{
		cerr << "Invalid p-value: Must be between 0 and 1.\n";
		exit(1);
	}
    ++optind;

    std::ifstream ifs(argv[optind]);
    if(!ifs) {
    	cerr <<"Could not open file "<<argv[optind]<<"\n";
    	exit(1);
    }
    SeqStreamSource seq_source(ifs);
    if (seq_source.eof())
    {
        cerr << "Error reading sequence\n";
        exit(1);
    }
    SeqIterator seq_it(&seq_source, buffer_size);

    ++optind;



    if (flat_bg)
        bg = flatBG(4);
    else if(bg.empty())
        bg = bgFromSequence(seq_it, 4, 0.1);

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
        if (matrix.size() == 4)
        {
            matrices.push_back(  counts2LogOdds( matrix, bg, 0.1 )  );
            thresholds.push_back(tresholdFromP(matrices.back(), bg, tol));
            names.push_back(&(argv[i]));
            cerr << "loaded matrix " << argv[i] << "\n";
        }
        else {
            cerr << "Matrix " << argv[i] << " has wrong alphabet size. Omitting.\n";
        }
    }

    cerr << "loaded matrices\n";

    int num_relevant_matrices = 0;
    for (unsigned int i = 0; i < matrices.size(); ++i)
    {
    	if(maxScore(matrices[i]) >= thresholds[i]) {
    		num_relevant_matrices++;
    	}
    }
    cerr<<"Matrices with possible hits: "<<num_relevant_matrices<<"\n";
    seq_it.reset();
    total = clock();
    position_t total_results = 0;

    if(algorithm < 4)
    {
        for (unsigned int i = 0; i < matrices.size(); ++i)
        {
            cout << *(names[i]) << "\n";
            matchArray results;

            seq_it.reset();
            single = clock();
            switch (algorithm)
            {
                case 0: results = naiveAlgorithm(seq_it, matrices[i], thresholds[i]); break;
                case 1: results = naiveSuperalphabetAlgorithmDNA(q, seq_it, matrices[i], thresholds[i]); break;
                case 2: results = permutatedLookAhead(seq_it, matrices[i], bg, thresholds[i]); break;
                case 3: results = lookaheadFiltrationDNA(q, seq_it, matrices[i], bg, thresholds[i]); break;

            }

            single = clock() - single;

            cerr << "Time: " << std::setprecision(5) << (double)single/((double)CLOCKS_PER_SEC) << "\n";
            cerr << "Hits: " << results.size() << "\n";
            total_results += results.size();
            if(print_results)
            	printMatchArray(results);
            cout << "\n";
        }
    }
    else
    {
        std::vector<matchArray> results;
        switch (algorithm)
        {
            case 4: results = multipleMatrixLookaheadFiltrationDNA(q, seq_it, matrices, bg, thresholds); break;
            case 5: results = searchDNA(seq_it, matrices, thresholds, bg, q); break;
        }

		for (unsigned int i = 0; i < matrices.size(); ++i)
		{
			cerr << *(names[i]) << "\n";
			cerr << "Hits: " << results[i].size() << "\n";
			total_results += results[i].size();
			if(print_results)
				printMatchArray(results[i]);
			cout << "\n";
		}
    }

    total = clock() - total;

    cerr << "Total time: " << std::setprecision(5) << (double)total/((double)CLOCKS_PER_SEC) << "\n";
    cerr << "Total results: " << total_results <<"\n";
    cout << "\n";


    return 0;
}


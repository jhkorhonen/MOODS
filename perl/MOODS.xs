// MOODS.xs - algorithm for finding PWM matches from sequences
// Copyright (C) 2009 Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.

#include <pssm_algorithms.hpp>
#include "seq_buffer.h"
#include <sstream>

using std::vector;

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"


/**
 * Calls perl object method length
 */
int call_length(SV *ref) {
	int i;
	dSP;
	ENTER;
	SAVETMPS;
	PUSHMARK(SP);
	XPUSHs(ref);
	PUTBACK;
	call_method("length", G_SCALAR);
	SPAGAIN;
	i = POPi;
	PUTBACK;
	FREETMPS;
	LEAVE;
	return i;
}

enum alphabet_type {
	ALPHABET_DNA,
	ALPHABET_RNA,
	ALPHABET_PROTEIN,
	ALPHABET_UNKNOWN
};

/**
 * Calls bio::seq method alphabet
 */
alphabet_type call_alphabet(SV *ref) {
	alphabet_type type;
	SV *alp;
	const char *s;
	STRLEN len;
	dSP;
	ENTER;
	SAVETMPS;
	PUSHMARK(SP);
	XPUSHs(ref);
	PUTBACK;
	call_method("alphabet", G_SCALAR);
	SPAGAIN;
	alp = POPs;
	s = SvPV(alp, len);
	if(strcmp(s, "dna")==0)
		type = ALPHABET_DNA;
	else if(strcmp(s, "rna")==0)
		type = ALPHABET_RNA;
	else if(strcmp(s, "protein")==0)
		type = ALPHABET_PROTEIN;
	else
		type = ALPHABET_UNKNOWN;

	PUTBACK;
	FREETMPS;
	LEAVE;
	return type;
}

/**
 * Converts perl double array to c++ equivalent
 */
doubleArray atoDoubleArraySV(SV *sv_array) {
	int i;
	int sv_size;
	SV **value;
	doubleArray c_array;

	if(SvTYPE(SvRV(sv_array)) != SVt_PVAV)
        return doubleArray();

    if(((sv_size = av_len((AV *) SvRV(sv_array))) < 0)) {
    	return doubleArray();
    }
    for(i = 0; i <= sv_size; i++) {
    	value = av_fetch((AV*) SvRV(sv_array), i, 1);
    	c_array.push_back(SvNV(*value));
    }
    return c_array;
}

/**
 * Converts a perl integer array reference to c++ object
 */
intArray atoIntArraySV(SV *sv_array) {
	int i;
	int sv_size;
	SV **value;
	intArray c_array;
	if(SvTYPE(SvRV(sv_array)) != SVt_PVAV)
        return intArray();

    if(((sv_size = av_len((AV *) SvRV(sv_array))) < 0)) {
    	return intArray();
    }
    for(i = 0; i <= sv_size; i++) {
    	value = av_fetch((AV*) SvRV(sv_array), i, 1);
    	c_array.push_back(SvIV(*value));
    }
    return c_array;
}

scoreArray atoScoreArraySV(SV *sv_array) {
	int i;
	int sv_size;
	SV **value;
	scoreArray c_array;

	if(SvTYPE(SvRV(sv_array)) != SVt_PVAV)
        return scoreArray();

    if(((sv_size = av_len((AV *) SvRV(sv_array))) < 0)) {
    	return scoreArray();
    }
    for(i = 0; i <= sv_size; i++) {
    	value = av_fetch((AV*) SvRV(sv_array), i, 1);
    	c_array.push_back(SvNV(*value));
    }
    return c_array;
}

scoreMatrix atoScoreMatrixSV(SV *sv_matrix) {
	int i, j, n, m;
	SV *row;
	SV **value;
	scoreMatrix c_matrix;

	if(SvTYPE(SvRV(sv_matrix)) != SVt_PVAV)
        return scoreMatrix();
	//Get number of rows
	if((n = av_len((AV *) SvRV(sv_matrix))) <= 0) {
		return scoreMatrix();
	}

	for(i = 0; i <= n; i++) {
        row = *av_fetch((AV *) SvRV(sv_matrix), i, 1);
        c_matrix.push_back(atoScoreArraySV(row));
    }
    return c_matrix;
}

class SeqSVSource : public SeqSourceI {
	char c;
	int value;
	SV *seq_ref;
	unsigned long offset;
	unsigned long seq_size;
public:

	SeqSVSource(SV *bio_seq) : seq_ref(bio_seq) {
		seq_size = call_length(seq_ref);
		offset = 0;
	}
	int read_data(char *dst, int maxlen) {
		const char *data = NULL;
		STRLEN length;
		SV *sv_string;
		int i;
		dSP;
		ENTER;
		SAVETMPS;
		unsigned long l_maxlen = maxlen;
		if((seq_size - offset) < l_maxlen) {
			l_maxlen = seq_size - offset;
			maxlen = (int) l_maxlen;
		}

		PUSHMARK(SP);
		XPUSHs(seq_ref);
		XPUSHs(sv_2mortal(newSViv(offset + 1)));
	    XPUSHs(sv_2mortal(newSViv(offset + maxlen)));
	    PUTBACK;
	    call_method("subseq", G_SCALAR);
		SPAGAIN;
		sv_string = POPs;
		data = SvPV(sv_string, length);
		int value = 0;
		alphabet_type type = call_alphabet(seq_ref);
		if(type == ALPHABET_DNA || type == ALPHABET_RNA) {
			for(i = 0; i < length; i++) {
				switch (data[i])
			    {
			        case 'a':
			        case 'A': value = 0; break;
			        case 'c':
			        case 'C': value = 1; break;
			        case 'g':
			        case 'G': value = 2; break;
			        case 'u':
			        case 'U':
			        case 't':
			        case 'T': value = 3; break;
			        case 'N':
			        default: value = (int) (4 * rand() / (1.0 + RAND_MAX)); break;
			    }
			    dst[i] = value;
			}
		}
		else {
			for(i = 0; i < length; i++) {
				value = data[i] - 'A';
				if (value >= 0 && value < 26)
					dst[i] = value;
			}
		}
		PUTBACK;
		FREETMPS;
	    LEAVE;
	    offset +=maxlen;
		return maxlen;
	}
	void reset() {
		offset = 0;
	}
	bool eof() {
		return offset >= seq_size;
	}
};

MODULE = MOODS		PACKAGE = MOODS

SV *
_search(seq, matrices, thresholds, algorithm, q, bgtype, combine, count_log_odds, threshold_from_p, buffer_size, sv_bg)
		SV *seq
        SV *matrices
        SV *thresholds
        int algorithm
        int q
        int bgtype
        bool combine
        bool count_log_odds
        bool threshold_from_p
        int buffer_size
        SV *sv_bg
    INIT:
    	int seq_length = call_length(seq);
    	alphabet_type type = call_alphabet(seq);
    	int alphabet_size = (type == ALPHABET_DNA || type == ALPHABET_RNA) ? 4 : 25;
   		SeqSVSource seq_source(seq);
   		SeqIterator seq_it(&seq_source, buffer_size >= 0 ? buffer_size : (seq_length > 100000000 ? 100000000 : seq_length));
        vector<scoreMatrix> c_matrices;
        vector<matchArray> c_ret;
        doubleArray d_thresholds;
        scoreArray c_thresholds;
        doubleArray bg;
        std::stringstream ss;
        int i;
        int num_matrices;
        SV *sv_matrix;
        AV *results;
        matchArray::iterator it;
        if (!SvROK(matrices) || !SvROK(thresholds)) {
            XSRETURN_UNDEF;
        }
    PPCODE:
    	if (bgtype == 3 && sv_bg) {
    		bg = atoDoubleArraySV(sv_bg);
    	}
    	else if (bgtype == 1) {
    	    bg = flatBG(alphabet_size);
    	}
    	else {
       		bg = bgFromSequence(seq_it, alphabet_size, 0.1);
    	}


    	num_matrices = av_len((AV *) SvRV(matrices)) + 1;

    	d_thresholds = atoDoubleArraySV(thresholds);

    	for(i = 0; i < num_matrices; i++) {
    		sv_matrix = *av_fetch((AV *) SvRV(matrices), i, 1);
    		if(count_log_odds)
    			c_matrices.push_back(counts2LogOdds(atoScoreMatrixSV(sv_matrix), bg, 0.1));
    		else
    			c_matrices.push_back(atoScoreMatrixSV(sv_matrix));

    		if(threshold_from_p)
    			c_thresholds.push_back(tresholdFromP(c_matrices.back(), bg, d_thresholds[i]));
    		else
    			c_thresholds.push_back(d_thresholds[i]);
    	}


       	if(type != ALPHABET_DNA && type != ALPHABET_RNA) {
       		for(i = 0; i < num_matrices; i++) {
       			seq_it.reset();
       			c_ret.push_back(naiveAlgorithm(seq_it, c_matrices[i], c_thresholds[i]));
       		}
       	}
        if(algorithm < 0 || algorithm > 5) {
       		c_ret = searchDNA(seq_it, c_matrices, c_thresholds, bg, q);
       	}
       	else if((num_matrices > 1) && (algorithm > 3) && combine) {
       		seq_it.reset();
       		switch(algorithm) {
       			case 4: //c_ret = multipleMatrixAhoCorasickLookaheadFiltrationDNA(q, c_seq, c_matrices, bg, c_thresholds); break;
       			case 5: c_ret = multipleMatrixLookaheadFiltrationDNA(q, seq_it, c_matrices, bg, c_thresholds); break;
       		};
       	}
       	else {
       		for(i = 0; i < num_matrices; i++) {
       			if(maxScore(c_matrices[i]) >= c_thresholds[i]) {
					seq_it.reset();
					switch(algorithm) {
						case 0: c_ret.push_back(naiveAlgorithm(seq_it, c_matrices[i], c_thresholds[i])); break;
						case 1: c_ret.push_back(naiveSuperalphabetAlgorithmDNA(q, seq_it, c_matrices[i], c_thresholds[i])); break;
						case 2: c_ret.push_back(permutatedLookAhead(seq_it, c_matrices[i], bg, c_thresholds[i])); break;
						case 3: //c_ret.push_back(ahoCorasickExpansionDNA(c_seq, c_matrices[i], c_thresholds[i])); break;
						case 4: //c_ret.push_back(ahoCorasickLookaheadFiltrationDNA(q, c_seq, c_matrices[i], bg, c_thresholds[i])); break;
						case 5: c_ret.push_back(lookaheadFiltrationDNA(q, seq_it, c_matrices[i], bg, c_thresholds[i])); break;
					};
       			}
       			else
       				c_ret.push_back(matchArray());
       		}
       	}

       	EXTEND(SP, c_ret.size());

       	for(i = 0; i < num_matrices; i++) {
       		results = (AV *)sv_2mortal((SV *)newAV());
       		for(it = c_ret[i].begin(); it != c_ret[i].end(); it++) {
            	av_push(results, newSVuv((*it).position));
            	av_push(results, newSVnv((*it).score));
        	}
       		PUSHs(newRV((SV *)results));
       	}


MODULE = MOODS		PACKAGE = MOODS

SV *
_count_log_odds(sv_matrix, sv_bg, ps)
		SV *sv_matrix
		SV *sv_bg
        double ps
    INIT:
        doubleArray bg;
        int i, k;
        scoreMatrix s_matrix = atoScoreMatrixSV(sv_matrix);
    PPCODE:
    	bg = atoDoubleArraySV(sv_bg);
   
  		scoreMatrix log_odds = counts2LogOdds(s_matrix, bg, ps);
  		
       	EXTEND(SP, log_odds.size());

		AV * results;
       	for(i = 0; i < log_odds.size(); i++) {
       		results = (AV *)sv_2mortal((SV *)newAV());
       		for(k = 0; k < log_odds[i].size(); k++) {
            	av_push(results, newSVnv(log_odds[i][k]));
        	}
       		PUSHs(newRV((SV *)results));
       	}
       	
       	
       	
MODULE = MOODS		PACKAGE = MOODS

SV *
_threshold_from_p(sv_matrix, sv_bg, p)
		SV *sv_matrix
		SV *sv_bg
        double p
    INIT:
        doubleArray bg;
        scoreMatrix matrix = atoScoreMatrixSV(sv_matrix);
    PPCODE:
    	bg = atoDoubleArraySV(sv_bg);
   
       	EXTEND(SP, 1);
       	PUSHs(newSVnv(tresholdFromP(matrix, bg, p)));


MODULE = MOODS		PACKAGE = MOODS

SV *
_bg_from_sequence(seq, ps)
		SV *seq
        double ps
    INIT:
        doubleArray bg;
        SeqSVSource seq_source(seq);
        SeqIterator seq_it(&seq_source, 1024);
        int i;
    PPCODE:
    	
       	bg = bgFromSequence(seq_it, 4, ps);

       	EXTEND(SP, bg.size());

       	for(i = 0; i <bg.size(); i++) {
       		PUSHs(newSVnv(bg[i]));
       	}

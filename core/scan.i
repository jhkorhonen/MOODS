%module scan

%{
#include "moods.h"
#include "moods_scan.h"
%}

%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(vector_double) vector<double>;
    %template(vector_uchar) vector<unsigned char>;
    %template(vector_vector_double) vector<vector<double> >;
    %template(vector_vector_vector_double) vector<vector<vector<double> > >;
    %template(vector_match) vector<MOODS::scan::match>;
    %template(vector_vector_match) vector<vector<MOODS::scan::match> >;
    %template(vector_string) vector<string>;
};

%include "moods.h"
%include "moods_scan.h"
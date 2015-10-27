%module parsers

%{
#include "moods.h"
#include "moods_parsers.h"
%}

%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(vector_double) vector<double>;
    %template(vector_uchar) vector<unsigned char>;
    %template(vector_vector_double) vector<vector<double> >;
};

%include "moods.h"
%include "moods_parsers.h"
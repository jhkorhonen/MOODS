
#ifdef SWIGPERL
%module "MOODS::tools"
#else
%module tools
#endif

%{
#include "moods.h"
#include "moods_tools.h"
#include "match_types.h"
%}

%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(vector_double) vector<double>;
    %template(vector_uchar) vector<unsigned char>;
    %template(vector_vector_double) vector<vector<double> >;
    %template(vector_variant) vector<MOODS::variant>;
};

%include "moods.h"
%include "match_types.h"
%include "moods_tools.h"
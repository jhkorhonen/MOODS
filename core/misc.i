
#ifdef SWIGPERL
%module "MOODS::misc"
#else
%module misc
#endif

%include "std_string.i"
%{
#include <cstdint>
#include "moods.h"
#include "moods_misc.h"
%}

%include "moods.h"
%include "moods_misc.h"
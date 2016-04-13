
#ifdef SWIGPERL
%module "MOODS::scan"
#else
%module scan
#endif


%{
#include "moods.h"
#include "moods_scan.h"
#include "match_types.h"
#include "scanner.h"
%}

%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(vector_double) vector<double>;
    %template(vector_uchar) vector<unsigned char>;
    %template(vector_size_t) vector<size_t>;
    %template(vector_vector_double) vector<vector<double> >;
    %template(vector_vector_vector_double) vector<vector<vector<double> > >;
    %template(vector_match) vector<MOODS::match>;
    %template(vector_vector_match) vector<vector<MOODS::match> >;
    %template(vector_string) vector<string>;
    %template(vector_variant) vector<MOODS::variant>;
    %template(vector_match_variant) vector<MOODS::match_with_variant>;
    %template(vector_vector_match_variant) vector<vector<MOODS::match_with_variant> >;
};

%include "moods.h"
%include "moods_scan.h"
%include "match_types.h"

class MOODS::scan::Scanner {
    public:
        Scanner(unsigned int window_size);
        Scanner(unsigned int window_size, const std::vector<std::string>& alphabet);

        void set_motifs(const std::vector<score_matrix>& matrices,
                        const std::vector<double>& bg,
                        const std::vector<double> thresholds);

        std::vector<std::vector<MOODS::match> > scan(const std::string& s);
        std::vector<std::vector<MOODS::match> > scan_max_hits(const std::string& s, size_t max_hits);

        std::vector<std::vector<MOODS::match_with_variant> > variant_matches(const std::string& seq, const std::vector<MOODS::variant> variants, int max_depth = 0);

};



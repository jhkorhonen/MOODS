%module scan

%{
#include "moods.h"
#include "moods_scan.h"
#include "scanner.h"
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

class MOODS::scan::Scanner {
    public:
        Scanner(unsigned int window_size);
        Scanner(unsigned int window_size, const std::vector<std::string>& alphabet);

        void set_motifs(const std::vector<MOODS::scan::Motif>& motifs);
        void set_motifs(const std::vector<score_matrix>& matrices,
                        const std::vector<double>& bg,
                        const std::vector<double> thresholds);

        std::vector<std::vector<scan::match> > scan(const std::string& s);

};



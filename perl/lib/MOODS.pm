# Module for MOODS

=head1 NAME

MOODS - Perl extension for finding significant matches of position weight matrices.

=cut

package MOODS;

use 5.010000;
use strict;
use warnings;
use Switch;
require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = qw(search);


our $VERSION = '0.01';

require XSLoader;
XSLoader::load('MOODS', $VERSION);

=head1 Subroutines

=head2 search

  Title   : search
  Usage   : @results = MOODS::search(seq => Bio::Seq(..), -matrix =>[[1,0],[0,1]] -threshold => 0.1)
  Function: Finds position weight matrix matches in dna sequence. 
  Returns : An array of references to result arrays. There are one result array
            corresponding to each matrix. (matrix1_results, matrix2_results,..)
            Each result array is a list of positions and scores like: 
            (pos1, score1, pos2, score2 ...) 
  Args    : 
           Obligatory
            -seq  BioPerl sequence object
            -matrix or -matrices
                 A matrix or a list of matrices. One matrix is represented 
                 as a typical perl multidimensional array: a reference to array of 
                 references to arrays of numbers, corresponding to the frequencies
                 or scores of the nucleotides A, C, G and T, respectively
            -threshold or -thresholds
                 A number or a list of numbers, used as threshold values for
                 matrix scanning.  If a single number is given, it is used
                 for all matrices; otherwise, there should be as many
                 threshold values as there are matrices.

           Optional
            -bg  Background distribution - an array of four doubles. If neither
                 -bg or -flatbg is given, the background is estimated from
                 the sequence. 
            -flatbg
                 If 1, the background distribution is set to a distribution
                 giving equal probability to all characters. Not compatible
                 with -bg. If neither -bg or -flatbg is given, the background
                 is estimated from the sequence. 
            -count_log_odds
                 If 1, assumes that the input matrices are frequency or
                 count matrices, and converts them to log-odds scoring
                 matrices; otherwise, treat them as scoring matrices.
                 Default 1.
            -threshold_from_p
                 If 1, assumes that thresholds are p-values and computes
                 the corresponding absolute threshold based on the matrix;
                 otherwise the threshold is used as a hard cut-off.
                 Default 1.
            -log_base
                 Base for logarithms used in log-odds computations. Relevant
                 if using -convert_log_odds => 1 and -threshold_from_p => 0.
                 Defaults to natural logarithm if parameter is not given.
            -pseudocount
                 Pseudocount used in log-odds conversion and added to
                 sequence symbol counts when estimating the background
                 from sequence. Default 1.

           Tuning parameters:
            (Optional, do not affect the results, but can give minor
             speed-ups in some cases. You can pretty much ignore these.)
            -algorithm  Selects the algorithm to use for scanning
                 "naive" naive algorithm
                 "pla" permutated lookahead algorithm
                 "supera" super alphabet algorithm. 
                   - Good for long matrices (> 20)
                 "lf" lookahead filtration algorithm. 
                   - Default algorithm in most cases.
                   - Sequence can be searched with multiple matrices 
                     simultaneously. 
                   - You should use this when you have large amount of matrices.
            -q   An integer, used for fine-tuning "supera" and "lf" algorithms.
                 The default value 7 should be ok pretty much always, but can 
                 be tuned to possibly slightly increase performance. 
            -combine
                 determines whether "lf" algorithm combines all
                 matrices to a single scanning pass. 
            -buffer_size

=cut

sub search {
    my %args = @_;
    
    #handle compulsory parameters
    if (not exists $args{-seq}) { die ("Sequence not specified\n"); }
    if(exists $args{-matrices} && (exists $args{-thresholds} || $args{-threshold})) {
   		if((ref $args{-matrices} ne "ARRAY"))  { die("Matrices argument must be an array"); }
   		for my $matrix (@{$args{-matrices}}) {
   			if((ref $matrix ne "ARRAY") || scalar(@{$matrix}) <= 0 ) { die("Invalid matrix parameter"); }
   		}
   		if(not exists $args{-thresholds}) {
   			$args{-thresholds} = [];
   			for (1..scalar(@{$args{-matrices}})) {
   				push(@{$args{-thresholds}}, $args{-threshold});
   			}
   			delete $args{-threshold};
   		}
   		elsif(scalar(@{$args{-matrices}}) != scalar(@{$args{-thresholds}})) {
   			die("Invalid parameters -matrices and/or -thresholds");
   		}
   	}
   	elsif (exists $args{-matrix} && exists $args{-threshold}) {
   		if((ref $args{-matrix} ne "ARRAY") || scalar(@{$args{-matrix}}) <= 0 ) {
   			die("Invalid matrix parameter");
   		}
   		$args{-matrices} = [$args{-matrix}];
   		$args{-thresholds} = [$args{-threshold}];
   		delete $args{-matrix};
   		delete $args{-threshold};
   	}
   	else {
   		die("Either -matrix or -matrices must be specified");
   	}
    
    # if(exists $args{-absolute_threshold} && $args{-absolute_threshold}) {
    #     $args{-count_log_odds} = 0;
    #     $args{-threshold_from_p} = 0;
    #     delete $args{-absolute_threshold};
    # }
    
    if(exists $args{-bg} && $args{-bg}) {
    	if((ref $args{-bg} ne "ARRAY") || scalar(@{$args{-bg}}) != 4)  { die("Invalid background"); }
    	$args{-bgtype} = 3;
    	if(exists $args{-flatbg}) { die("You can't specify both bg and flatbg params"); }
    }
    elsif(exists $args{-flatbg} && $args{-flatbg}) {
    	$args{-bgtype} = 1;
    	delete $args{-flatbg};
    }
    
    #default parameters
    my %params = (-seq => 0, -matrices => 0, -thresholds => 0, -algorithm => "undefined", -q => 7, -bgtype => 0, -combine => 1, -count_log_odds => 1, -count_log_odds => 1, -threshold_from_p => 1, -buffer_size => -1, -bg => 0, -pseudocount => 1, -log_base => 0);
    
    #check for invalid parameters
    foreach my $param (keys %args) {
    	if(not exists $params{$param}) { die("Invalid parameter: ".$param); }
    	$params{$param} = $args{$param};
    }
   
  
   	my $alg = 4;
   	switch($params{-algorithm}) {
   		case "naive" { $alg = 0 }
   		case "supera" { $alg = 1 }
   		case "pla" { $alg = 2 }
   		#case "ace" { $alg = 3 }
   		#case "acelf" { $alg = 4 }
   		case "lf" { $alg = 5 }
   		case "undefined" { $alg = -1 }
   		else { die("Invalid -algorithm parameter"); }
   	}
   	
   	my @ret = _search($params{-seq}, $params{-matrices}, $params{-thresholds},$alg, $params{-q}, $params{-bgtype}, $params{-combine}, $params{-count_log_odds}, $params{-threshold_from_p}, $params{-buffer_size}, $params{-bg}, $params{-pseudocount}, $params{-log_base});
   	
   	return wantarray ? @ret : \@ret;
}


sub countLogOdds {
	my @ret = _count_log_odds(@_);
	return wantarray ? @ret : \@ret;
}

sub thresholdFromP {
	return _threshold_from_p(@_);
}

sub bgFromSequence {
	my @ret = _bg_from_sequence(@_);
	return wantarray ? @ret : \@ret;
}

1;

__END__

=head1 SYNOPSIS

  use Bio::Perl;
  use Bio::Seq;
  use MOODS;
  use MOODS::Tools qw(printResults);
  
  #we need a position weight matrix
  my $matrix = [ [10,0,0],
                 [0,10,0],
                 [0,0,10],
                 [10,10,10]];
  
  #we need also a bioperl sequence object
  my $seq = Bio::Seq->new(-seq              => 'actgtggggacgtcagtagcaggcatag',
                          -alphabet         => 'dna' );
                          
  my @results = MOODS::search(-seq => $seq, -matrix => $matrix, -threshold => 0.3);
  
  printResults($results[0]);

=head1 SEE ALSO

BioPerl documentation.

=head1 AUTHOR

Petri J Martinmaki, Janne H Korhonen

=head1 COPYRIGHT AND LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/

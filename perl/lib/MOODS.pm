# Module for MOODS
#
# Copyright Petri Martinmäki

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
               You can give one or multiple matrices. One matrix is represented 
               as a typical perl multidimensional array: a reference to array of 
               references to arrays of numbers.
            -threshold or -thresholds
               One threshold is ok for as many matrices as you want, but you
               can give different threshold for each matrix. 
           Optional
            -algorithm  You can switch search algorithm (Doesn't change results)
                "naive" naive algorithm
                "pla" permutated lookahead algorithm
                "supera" super alphabet algorithm. 
                  - Good for long matrices (> 20)
                "lf" lookahead filtration algorithm. 
                  - Default algorithm in most cases.
                  - Sequence can be searched with multiple matrices 
                    simultaneously. 
                  - You should use this when you have large amount of matrices.
            -q  You can optionally give a parameter for algorithm (Doesn't 
                change results)
            -combine  If you use lookahead filtration algorithm and don't wan't
                to use multiple matrix simultaneous search, you can spesify
                this parameter to 0. (usually slows down)
            -absolute_threshold  1 if threshold is given as an absolute value.
            -buffer_size
            -bg  Background distribution - an array of four doubles. By default
                 the background is calclulated from sequence. 
            -flatbg  1 if background distribution is not calculated from 
                 the sequence.

=cut

sub search {
    my %args = @_;
    
    #handle compulsory parameters
    if (not exists $args{-seq}) { die ("Sequence not spesified\n"); }
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
   		die("Either -matrix or -matrices must be spesified");
   	}
    
    if(exists $args{-absolute_threshold} && $args{-absolute_threshold}) {
    	$args{-count_log_odds} = 0;
    	$args{-threshold_from_p} = 0;
    	delete $args{-absolute_threshold};
    }
    
    if(exists $args{-bg} && $args{-bg}) {
    	if((ref $args{-bg} ne "ARRAY") || scalar(@{$args{-bg}}) != 4)  { die("Invalid background"); }
    	$args{-bgtype} = 3;
    	if(exists $args{-flatbg}) { die("You can't spesify both bg and flatbg params"); }
    }
    elsif(exists $args{-flatbg} && $args{-flatbg}) {
    	$args{-bgtype} = 1;
    	delete $args{-flatbg};
    }
    
    #default parameters
    my %params = (-seq => 0, -matrices => 0, -thresholds => 0, -algorithm => "undefined", -q => 7, -bgtype => 0, -combine => 1, -count_log_odds => 1, -threshold_from_p => 1, -buffer_size => -1, -bg => 0);
    
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
   	
   	my @ret = _search($params{-seq}, $params{-matrices}, $params{-thresholds},$alg, $params{-q}, $params{-bgtype}, $params{-combine}, $params{-count_log_odds}, $params{-threshold_from_p}, $params{-buffer_size}, $params{-bg});
   	
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

Petri J Martinmaki, E<lt>pmartinm@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Petri J Martinmaki

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

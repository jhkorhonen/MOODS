#!/usr/bin/perl

use warnings;
use strict;

use Bio::Perl;
use Bio::Seq;
use MOODS;
use MOODS::Tools qw(numResults);

#sequence loaded form file
my $seq_io = Bio::SeqIO->new(-file => "data/sequence/dnaACGT.txt", -format => "fasta");
my $seq = $seq_io->next_seq;


#we specify two matrices
my $matrix1 =
			[	[0,1,0,0,0,0,0,1,1,0],
				[1,0,0,0,0,0,0,0,0,0],
				[0,0,0,0,0,0,0,0,0,0],
				[0,0,1,1,1,1,1,0,0,1]
			];
my $matrix2 = 
			[	[10,0,10,3,5,5,0],
				[ 0,5, 0,3,5,0,5],
				[ 0,1, 0,3,0,5,0],
				[ 0,4, 0,1,0,0,5]
			];

# By default the matrices are converted to log-odds scoring matrices and the actual score threshold
# is computed from p-value given as -threshold 

my @results = MOODS::search(-seq => $seq, -matrices => [$matrix1, $matrix2], -thresholds => [0.001, 0.01], -flatbg => 1);


#we print only the number of results at this time
print "Matrix1 results: ".numResults($results[0])."\n";
print "Matrix2 results: ".numResults($results[1])."\n";
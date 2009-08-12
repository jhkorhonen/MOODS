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


#we spesify two matrices
my $matrix1 =
			[	[0,1,0,0,0,0,0,1,1,0],
				[1,0,0,0,0,0,0,0,0,0],
				[0,0,0,0,0,0,0,0,0,0],
				[0,0,1,1,1,1,1,0,0,1]
			];
my $matrix2 = 
			[	[10,0,10,3,5,5],
				[0,5,0,3,5,0,5],
				[0,1,0,3,0,5,0],
				[0,4,0,1,0,0,5]
			];

#we use the same threshold for both matrices, we could also spesify -thresholds => [0.011, 0.011]
my @results = MOODS::search(-seq => $seq, -matrices => [$matrix1, $matrix2],-threshold => 0.011, -flatbg => 1);


#we print only the number of results at this time
print "Matrix1 results: ".numResults($results[0])."\n";
print "Matrix2 results: ".numResults($results[1])."\n";
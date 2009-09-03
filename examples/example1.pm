#!/usr/bin/perl

use warnings;
use strict;

use Bio::Perl;
use Bio::Seq;
use MOODS qw(search);
use MOODS::Tools qw(printResults);

#we need a position weight matrix
my $matrix = [ [10,0,0],
            [0,10,0],
            [0,0,10],
            [10,10,10]];

#we need also a bioperl sequence object
my $seq = Bio::Seq->new(-seq              => 'actgtggcgtcaacgtaggccaacgtggacccgtacgtaaacgaagaggggtagtc',
                       -alphabet         => 'dna' );


my @results = search(-seq => $seq, -matrix => $matrix, -threshold => 30, -absolute => 1);

printResults($results[0]);
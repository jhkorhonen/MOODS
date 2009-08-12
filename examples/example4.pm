#!/usr/bin/perl

use warnings;
use strict;

use Bio::Perl;
use Bio::Seq;
use MOODS;
use MOODS::Tools qw(printResults numResults readMatrix);


my $seq_io = Bio::SeqIO->new(-file => "data/sequence/dnaACGT.txt", -format => "largefasta");
my $seq = $seq_io->next_seq;
print "loaded sequence\n";

my $matrix = readMatrix("data/matrix/JASPAR_CORE_2008/MA0110.pfm");
print "matrix loaded\n";

my @results = MOODS::search(-seq => $seq, -matrix => $matrix, -threshold => 0.001, -buffer_size => 600000);

#we only print first twenty results
printResults($results[0], [0..20]);

print "Sequence length: ".$seq->length."\n";
print "Total results: ".numResults($results[0])."\n";
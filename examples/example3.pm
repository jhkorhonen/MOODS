#!/usr/bin/perl

use warnings;
use strict;

use Bio::Perl;
use Bio::Seq;
use MOODS;
use MOODS::Tools qw(printResults numResults readMatrix);


my $seq_io = Bio::SeqIO->new(-file => "data/sequence/dnaACGT.txt", -format => "fasta");
my $seq = $seq_io->next_seq;
print "loaded sequence\n";

my $matrix = readMatrix("data/matrix/JASPAR_CORE_2008/MA0110.pfm");
print "matrix loaded\n";


# We convert the count matrix into log-odds scores in base 2 logarithm and find hits scoring
# better than 8 with this scoring matrix

my @results = MOODS::search(-seq => $seq, -matrix => $matrix, -threshold => 8, -threshold_from_p => 0,
                            -count_log_odds => 1, -log_base => 2);

printResults($results[0]);

print "Sequence length: ".$seq->length."\n";
print "Total results: ".numResults($results[0])."\n";
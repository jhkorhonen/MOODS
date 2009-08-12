#!/usr/bin/perl

use warnings;
use strict;

use Bio::Perl;
use Bio::Seq;
use MOODS;
use MOODS::Tools qw(numResults readMatrix);


#very simple parsing of command line parameters
my ($window_length, $threshold, $k_value, $sequence_file, @matrix_files) = @ARGV;

if((not defined($window_length)) || (not defined($threshold)) || (not defined($sequence_file)) || (not @matrix_files)) {
	print "Usage: <script> window_width threshold hit_threshold sequence matrix1 matrix2 ...\n";
	print "Example: \n";
	print "perl window_search.pm 100 0.005 0 data/sequence/dnaACGT.txt data/matrix/JASPAR_CORE_2008/MA0130.pfm data/matrix/JASPAR_CORE_2008/MA0124.pfm data/matrix/JASPAR_CORE_2008/MA0015.pfm\n";
	die("Invalid parameters");
}
#load a sequence from file
my $seq_io = Bio::SeqIO->new(-file => $sequence_file, -format => "fasta", -alphabet => "dna");
my $seq = $seq_io->next_seq;
my $seq_length = $seq->length;
print(">Sequence loaded.\n");


#load matrices from files
my @matrices;
foreach my $matrix_file (@matrix_files) {
	my $matrix = readMatrix($matrix_file);
	push(@matrices, $matrix);
}
print(">Matrices loaded.\n");


#actual search of PWM matches
my @results = MOODS::search(-seq => $seq, -matrices => [@matrices], -threshold => $threshold);
print(">Sequence searched\n");

my $num_matrices = scalar(@matrices);
my $last_matrix = $num_matrices - 1;
my @num_hits;
my $sum_hits;
my @matrix_sizes;
my @next_hit_indices;
my @last_hit_indices;
my $top_window_begin = 0;
my $top_window_end = 0;
my $top_hits = 0;
my $current_window_begin = 0;
my $current_window_end = 0;
my $current_hit_matrix = 0;
my $current_hits = 0;
my $last_hits = 0;

foreach my $i (0..$last_matrix) {
	push(@next_hit_indices, 0);
	push(@last_hit_indices, 0);
	push(@matrix_sizes, scalar(@{$matrices[$i]->[0]}));
	push(@num_hits, scalar(@{$results[$i]}) / 2);
	$sum_hits += scalar(@{$results[$i]}) / 2;
}
print "Sequence length: ".$seq_length."\n";
print "Total number of hits: ".$sum_hits."\n";
if($k_value) {
	print "Searching windows containing at least ".$k_value." hits:\n";
}


#loop through window positions
while(1) {
	my $i = 0;
	#move back of the window forward
	my $smallest = -1;	#smallest hit position not jet checked
	foreach $i (0..$last_matrix) {
		if($last_hit_indices[$i] < $num_hits[$i] && ($smallest == -1 || ($results[$i]->[$last_hit_indices[$i] * 2] < $smallest))) {
			$smallest = $results[$i]->[$last_hit_indices[$i] * 2];
		}
	}
	if($smallest == -1) {
		last; #all positions checked
	}
	$current_window_begin = $smallest;
	
	
	#add all hits that are inside the window
	$last_hits = $current_hits;
	foreach $i (0..$last_matrix) {
		while($next_hit_indices[$i] < $num_hits[$i] && $results[$i]->[$next_hit_indices[$i] * 2] - $current_window_begin + $matrix_sizes[$i] < $window_length) {
			if($current_window_end < $results[$i]->[$next_hit_indices[$i] * 2] + $matrix_sizes[$i]) {
 				$current_window_end = $results[$i]->[$next_hit_indices[$i] * 2] + $matrix_sizes[$i];
			}
			$current_hits++;
			if($k_value && $current_hits >= $k_value) {
				print "Hits ".$current_hits." matrix size: ".$matrix_sizes[$i]." hit position: ".$results[$i]->[$next_hit_indices[$i] * 2]."\n";
			}
			$next_hit_indices[$i]++;
		}
	}
	
	#check if current is the best window
	if($current_hits >= $top_hits) {
		$top_hits = $current_hits;
		$top_window_begin = $current_window_begin;
		$top_window_end = $current_window_end;
	}
	
	#print all windows where number of hits is greater than k_value
	if($k_value && $current_hits >= $k_value && $current_hits > $last_hits) {
		print "\tlength: ".($current_window_end - $current_window_begin)."\n";
		print "\tposition: [".$current_window_begin.", ".($current_window_end-1)."]\n";
		print "\thits: ".$current_hits."\n";
	}
	
	
	#remove hits, that are not inside the window on next round
	foreach $i (0..$last_matrix) {
		if($last_hit_indices[$i] < $num_hits[$i] && $results[$i]->[$last_hit_indices[$i] * 2] == $smallest) {
			$last_hit_indices[$i]++;
			$current_hits--;
		}
	}
}
print "Window containing most hits: \n";
if($top_hits) {
	print "\tlength: ".($top_window_end - $top_window_begin)."\n";
	print "\tposition: [".$top_window_begin.", ".($top_window_end-1)."]\n";
	print "\thits: ".$top_hits."\n";
}

1;
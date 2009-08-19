#!/usr/bin/perl

use warnings;
use strict;



my $list ="";
my ($manifest_filename) = @ARGV;
open(FILE, "<", $manifest_filename);
my @folders = ("");
my $num_white = 0;
for my $line (<FILE>) {
	chop($line);
	$num_white = 0;
	while ($line =~ s/ (.+)/$1/) {
	    $num_white++;
	}
	if($num_white %4 != 0) { die("Invalid indention at line: ".$line); }
	$num_white = $num_white/4;
	$folders[$num_white + 1] = $folders[$num_white].$line;
	my $filename = "../../MOODS/$folders[$num_white + 1]";
	-e $filename or die ("file not found: ".$filename);
	if(not $filename =~ m/.+\/$/) {
		$list .= " ".$filename;	
		print $filename."\n";
	}
}
close(FILE);

system("tar cvvzf MOODS-v1.0.tar.gz ".$list) == 0 or die("creating package failed");

1;
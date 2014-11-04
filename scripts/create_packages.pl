#!/usr/bin/perl

use warnings;
use strict;

sub create_package {
	my $manifest_filename = shift;
	my $package_name = shift;
	my $list ="";
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
		my $filename = "MOODS/$folders[$num_white + 1]";
		-e "../../$filename" or die ("file not found: ".$filename);
		if(not $filename =~ m/.+\/$/) {
			$list .= " ".$filename;	
			print $filename."\n";
		}
	}
	close(FILE);
	system("COPYFILE_DISABLE=1 tar --exclude *.DS_Store -C ../../ -cvvzf ".$package_name." ".$list) == 0 or die("creating package failed");
}

create_package("package.txt", "MOODS.tar.gz");

1;
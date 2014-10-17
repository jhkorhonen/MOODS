#!/usr/bin/perl

use warnings;
use strict;

use Pod::Simple::HTML;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
 
sub convert {
	my $result;
	my $in_file = shift;
	my $out_file = shift;
	open OUT_FILE, ">".$out_file or die("Cannot open ".$out_file);
	my $parser = Pod::Simple::HTML->new();
	$parser->index(1);
	$parser->html_css('http://search.cpan.org/s/style.css');
	$parser->output_fh(*OUT_FILE);
	$parser->parse_file($in_file);
}

rmdir("doc");
mkdir("doc");
mkdir("doc/perl");

convert("../perl/lib/MOODS.pm", "doc/perl/moods.html");
convert("../perl/lib/MOODS/Tools.pm", "doc/perl/tools.html");


mkdir("doc/python");
rcopy("../python/MOODS/__init__.py", "doc/python/MOODS/") or die("../python/MOODS not found");
chdir("doc/python");
system("pydoc -w MOODS")== 0 or die("Generating python documents failed");

1;
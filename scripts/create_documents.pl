#!/usr/bin/perl

use warnings;
use strict;

use Pod::Simple::HTML;
 
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


convert("../perl/lib/MOODS.pm", "doc/moods.html");
convert("../perl/lib/MOODS/Tools.pm", "doc/tools.html");

1;
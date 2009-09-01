# Tools module for MOODS
#
# Copyright Petri Martinmäki

package MOODS::Tools;
use warnings;
use strict;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(readMatrix printMatrix numResults printResults getResult);


=head1 NAME

MOODS::Tools - Some useful tools to work with parameters and results of PWM search.

=head1 Subroutines

=cut

=head2 readMatrix

  Title   : readMatrix
  Usage   : my $matrix = readMatrix("matrixfile");
  Function: reads a matrix from file
  Returns : matrix
  Args    : filename

=cut

sub readMatrix {
	my ($filename) = @_;
	my $file;
	open($file, "<", $filename) or die "File not found".$filename;
	my @ret = ();
	my $line;
	for $line (<$file>) {
		chomp($line);
		$line =~ s/^\s+//;
		my @numbers = split(/[\s]+/, $line);
		push(@ret, \@numbers);
	} 
	return \@ret;
}

=head2 printMatrix

  Title   : printMatrix
  Usage   : my $matrix = readMatrix("matrixfile");
  			printMatrix($matrix)
  Function: prints matrix
  Args    : matrix reference

=cut

sub printMatrix {
	my $matrixref = shift;
	foreach my $row (@{$matrixref}) {
		foreach my $num (@{$row}) {
			print $num." ";
		}
		print "\n";
	}
}

=head2 reverseComplement

  Title   : reverseComplement
  Usage   : my $matrix2 = reverseComplement($matrix1);
  Args    : matrix reference

=cut

sub reverseComplement {
	my $source = shift;
	my @dest = ();
	for my $i (3,2,1,0) {
		my @tmp = reverse(@{$source->[$i]});
		push(@dest, \@tmp);
	}
	return \@dest;
}

=head2 numResults

  Title   : numResults
  Usage   : my @results = MOODS::search(..);
  			$num_results = numResults(@results[0]);
  Function: calculates a number of matches
  Args    : an array reference to result array

=cut

sub numResults {
	my $resultref = shift;
	my $num = @{$resultref};
	return $num/2;
}

=head2 getResult

  Title   : getResult
  Usage   : my @results = MOODS::search(..);
  			$num_results = numResults($results[0]);
  			foreach $i (0..$num_results) {
  				my ($pos, $score) = getResult($results[0], $i); 
  			}
  Function: gets a position score pair from a result array
  Args    : a result array and an index

=cut

sub getResult {
	my $resultref = shift;
	my $position = shift;
	return @{$resultref}[$position*2, $position*2+1];
}

=head2 printResults

  Title   : printResults
  Usage   : my @results = MOODS::search(..);
  			printResults($results[0]);
  Function: print a result array
  Args    : an array reference to result array

=cut

sub printResults {
	my $resultref = shift;
	my $positionref = shift;
	if((not defined($positionref)) || (ref $positionref ne "ARRAY")) {
		$positionref = [0..(numResults($resultref)-1)];
	}
	foreach my $i (@{$positionref}) {
		my ($position, $score) = getResult($resultref, $i);
		print "Position: ".$position." Score: ".$score."\n";
	}
}
1;

__END__

=head1 AUTHOR

Petri J Martinmaki, E<lt>pmartinm@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Petri J Martinmaki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/
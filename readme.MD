MOODS version 1.0.2
=======================

MOODS is a suite of algorithms for matching position weight matrices (PWM) against DNA sequences. It features advanced matrix matching algorithms implemented in C++ that can be used to scan hundreds of matrices against chromosome-sized sequences in few seconds.

MOODS can be used as standalone analysis tool or as a component in larger programs. It contains Python and Perl/BioPerl interfaces, and can thus be easily called from C++, Python and Perl programs.

MOODS is dual-licenced under GPL version 3 license and Biopython license. See COPYING.GPLv3 and COPYING.BIOPYTHON files for license details.

The project web page is at http://www.cs.helsinki.fi/group/pssmfind. Please refer there for further information and documentation.

Contact
-------

MOODS has been written by Pasi Rastas, Janne H Korhonen and Petri Martinmäki. The project is currently maintained by Janne H Korhonen (jazkorho (at) cs.helsinki.fi).


Installation
------------

You need to compile the C++ library before installing the Perl and Python interfaces. This can be done by running "make" in the "src" directory.

Installing extensions:

 * For perl extension see perl/README.
 * For python extension see python/README

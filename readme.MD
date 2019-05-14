MOODS 1.9.4
===========

MOODS is a suite of algorithms for matching position weight matrices (PWM) against DNA sequences, featuring advanced matrix matching algorithms implemented in C++ that can be used to scan hundreds of matrices against chromosome-sized sequences in few seconds. MOODS has been designed with integration into large-scale python workflows in mind, but can also be used as a stand-alone analysis tool.

 * MOODS 1.9.4 released:
   * PyPi package (https://pypi.org/project/MOODS-python/), courtesy of [Fabio Ticconi](https://github.com/fabio-t)
   * support for adjusting the precision of p-value computation
   * bugfixes for variant matching code
 * MOODS 1.9.3 released:
   * fixed bug in background distribution estimation with non-ACGT characters
   * updated install script for better support of newer macOS releases
 * new MOODS paper: [Fast motif matching revisited: high-order PWMs, SNPs and indels](https://doi.org/10.1093/bioinformatics/btw683)
 * MOODS 1.9.2 released:
   * fixed installation script
 * MOODS 1.9.1 released:
   * introduced python commandline tool for basic analysis
   * bugfixes in p-value and scanner code
   * minor changes in `MOODS.parsers` interface

The current release is the 1.9 series of MOODS, an intermediate release towards the fully-revamped 2.0 release later on. Currently all the algorithmic parts are fully ready for use, and should be preferred over the 1.0 releases when possible. However, the documentation and some input/ouput-related functions are still in development.

Download
--------

Download the package:

 * Python 2.7: [MOODS-python-1.9.4.tar.gz](https://github.com/jhkorhonen/MOODS/releases/download/v1.9.4/MOODS-python-1.9.4.tar.gz)
 
Or install from [PyPi](https://pypi.org/project/MOODS-python/) using pip:

 * `pip install --user MOODS-python`

Currently only the Python package is available; it contains MOODS Python libraries along a commandline interface for most of the MOODS functions.


Usage
-----

See the GitHub wiki for documentation:

 * [Installation](https://github.com/jhkorhonen/MOODS/wiki/Installation)
 * [Getting started](https://github.com/jhkorhonen/MOODS/wiki/Getting-started)


Contact
-------

 * For bug reports and other software issues, preferably use [issue tracker](https://github.com/jhkorhonen/MOODS/issues).
 * [Janne H. Korhonen](https://users.aalto.fi/~korhonj2/) MOODS developer
 
Licensing
---------

MOODS is dual-licenced under GPL version 3 license and Biopython license. See COPYING.GPLv3 and COPYING.BIOPYTHON files for license details.

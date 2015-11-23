MOODS 1.9.0
===========

MOODS is a suite of algorithms for matching position weight matrices (PWM) against DNA sequences. It features advanced matrix matching algorithms implemented in C++ that can be used to scan hundreds of matrices against chromosome-sized sequences in few seconds.

This is the 1.9 series of MOODS, an intermediate release towards the fully-revamped 2.0 release later on. Currently all the algorithmic parts are fully ready for use, and should be preferred over the 1.0 releases when possible. However, the documentation and some input/ouput-related functions are still in development.

Download
--------

Download the packages here:

* Python [MOODS-python-1.9.0.tar.gz](tarballs/MOODS-python-1.9.0.tar.gz)

Currently only the Python interfaces are available. Standalone command line application and Perl interfaces are in development.


Usage
-----

### Python

To start using the Python interfaces, untar the package and compile the C++ extensions:

    tar xzvf MOODS-python-1.9.0.tar.gz
    cd MOODS-python-1.9.0/
    python setup.py build_ext --inplace

The result should be a functional MOODS python package under the `MOODS/` subdirectory (tested on Linux with gcc 4.6.4 and on Max OS X with Clang 7.0.0; the C++ parts use C++11 features, so make sure your compiler is new enough).

You can test the package using the provided example scripts and example data; these also show how to use the various features of the package. For example, try:

    python ex-basic-usage.py example-data/seq/chr1-5k-55k.fa example-data/matrices/


Licensing
---------

MOODS is dual-licenced under GPL version 3 license and Biopython license. See COPYING.GPLv3 and COPYING.BIOPYTHON files for license details.

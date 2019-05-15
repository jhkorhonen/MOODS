MOODS-python 1.9.4.1
====================

Motif Occurrence Detection Suite

https://www.cs.helsinki.fi/group/pssmfind/


Installing from PyPi
--------------------

The quickest way to install MOODS is using pip:

(a) Install the module to your local Python library path:

    pip install --user MOODS-python

(b) Install the module for all users:

    sudo pip install MOODS-python

Installing the MOODS package from PyPi installs the `moods-dna.py` script to
your Python executable directory. Test it by running:

    moods-dna.py --help

(Note that this assumes that the appropriate Python executable directory is in
your path. For example, if you install the package with --user, you can do this
with

    PATH="$PATH:`python -m site --user-base`/bin"

if it doesn't work otherwise.)


Installing from Source
----------------------

Alternatively, you can install MOODS from the source package, for example if you
don't have your pip configured or you just want a quick throwaway install.

(a) If you don't want to bother installing MOODS, you can compile the extensions
    in place and create a symlink under `scripts/`:

    python setup.py build_ext --inplace
    cd scripts
    ln -s ../MOODS/

(b) Install the module to your local Python library path:

    python setup.py install --user

(c) Install the module for all users:

    sudo python setup.py install

You can test MOODS by running the `moods-dna.py` script with
the example data provided with the package:

    cd scripts/
    python moods-dna.py -m example-data/matrices/*.{pfm,adm} -s example-data/seq/chr1-5k-55k.fa -p 0.0001
    python moods-dna.py --help


Documentation
-------------

Documentation on GitHub wiki: https://github.com/jhkorhonen/MOODS/wiki


Copyright and Licence
---------------------

Copyright (C) 2009 by Petri J Martinmäki, Janne Korhonen, Pasi Rastas

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version, or under the terms of the Biopython
License.

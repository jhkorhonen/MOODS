MOODS-python 1.9.3
==================

Motif Occurrence Detection Suite

https://www.cs.helsinki.fi/group/pssmfind/


Getting Started
---------------

Start by compiling the MOODS extensions and optionally installing
the MOODS module to Python library path. Do one of the following:

(a) Compile the extensions in place and create a symlink under `scripts/`
(if you just want to test things without installing to Python library path,
or you are simply planning to use the `moods_dna.py` scripts for basic 
analysis):

    python setup.py build_ext --inplace
    cd scripts
    ln -s ../MOODS/    

(b) Install the module to your local Python library path:

    python setup.py install --user

(c) Install the module for all users:

    sudo python setup.py install

In all cases, you can test MOODS by running the `moods_dna.py` script
with the example data provided with the package:

    cd scripts/
    python moods_dna.py -m example-data/matrices/*.{pfm,adm} -s example-data/seq/chr1-5k-55k.fa -p 0.0001
    python moods_dna.py --help


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

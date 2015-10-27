#!/bin/bash
swig -c++ -python -outdir MOODS/ core/scan.i
swig -c++ -python -outdir MOODS/ core/tools.i
swig -c++ -python -outdir MOODS/ core/misc.i
swig -c++ -python -outdir MOODS/ core/parsers.i
python setup.py build_ext --inplace
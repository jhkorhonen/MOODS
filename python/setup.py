#!/usr/bin/env python

"""
setup.py file for MOODS
"""

from setuptools import setup, Extension
from os import path

common_includes = ["core/"]
common_compile_args = ['-march=native', '-O3', '-fPIC', '--std=c++11']

tools_mod = Extension('MOODS._tools',
                           sources=['core/tools_wrap.cxx',
                                    'core/moods_tools.cpp',
                                    'core/moods_misc.cpp',
                                    'core/match_types.cpp'],
                           include_dirs=common_includes,
                           extra_compile_args=common_compile_args,
                           )

scan_mod = Extension('MOODS._scan',
                           sources=['core/scan_wrap.cxx',
                                    'core/moods_scan.cpp',
                                    'core/motif_0.cpp',
                                    'core/motif_h.cpp',
                                    'core/moods_misc.cpp',
                                    'core/scanner.cpp',
                                    'core/moods_tools.cpp',
                                    'core/match_types.cpp'
                                ],
                           include_dirs=common_includes,
                           extra_compile_args=common_compile_args,
                           )

parsers_mod = Extension('MOODS._parsers',
                           sources=['core/parsers_wrap.cxx',
                                    'core/moods_parsers.cpp',
                                    'core/moods_misc.cpp',
                                    'core/moods_tools.cpp',
                                    'core/match_types.cpp'],
                           include_dirs=common_includes,
                           extra_compile_args=common_compile_args,
                           )



this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'readme.MD')) as f:
    long_description = f.read()

setup (name = 'MOODS-python',
       version = '1.9.4.1',
       description = 'MOODS: Motif Occurrence Detection Suite',
       long_description = long_description,
       long_description_content_type='text/markdown',
       maintainer = "Janne H. Korhonen",
       maintainer_email = "janne.h.korhonen@gmail.com",
       url='https://www.cs.helsinki.fi/group/pssmfind/',
       license = "GPLv3 / Biopython license",
       ext_modules = [tools_mod, scan_mod, parsers_mod],
       packages = ["MOODS"],
       scripts=['scripts/moods-dna.py'],
       classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
       keywords="PWM, PSSM, motif scan"
)

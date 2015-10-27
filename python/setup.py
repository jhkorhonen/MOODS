#!/usr/bin/env python

"""
setup.py file for MOODS
"""

from distutils.core import setup, Extension


tools_mod = Extension('MOODS._tools',
                           sources=['core/tools_wrap.cxx', 'core/moods_tools.cpp','core/moods_misc.cpp'],
                           include_dirs=["core/"],
                           extra_compile_args=['-march=native', '-O3', '-fPIC', '--std=c++0x'],
                           )

scan_mod = Extension('MOODS._scan',
                           sources=['core/scan_wrap.cxx',
                                    'core/moods_scan.cpp',
                                    'core/motif_0.cpp',
                                    'core/motif_h.cpp',
                                    'core/moods_misc.cpp',
                                    'core/scanner.cpp'],
                           include_dirs=["core/"],
                           extra_compile_args=['-march=native', '-O3', '-fPIC', '--std=c++0x'],
                           )

parsers_mod = Extension('MOODS._parsers',
                           sources=['core/parsers_wrap.cxx',
                                    'core/moods_parsers.cpp',
                                    'core/moods_misc.cpp',
                                    'core/moods_tools.cpp'],
                           include_dirs=["core/"],
                           extra_compile_args=['-march=native', '-O3', '-fPIC', '--std=c++0x'],
                           )

setup (name = 'MOODS',
       version = '1.9',
       description = 'MOODS: Motif Occurrence Detection Suite',
       ext_modules = [tools_mod, scan_mod, parsers_mod],
       # ext_modules = [tools_mod, misc_mod],
       py_modules = ["MOODS.tools", "MOODS.scan"],
)
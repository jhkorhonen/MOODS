       #!/usr/bin/env python

"""
setup.py file for haplotyping_core
"""

from distutils.core import setup, Extension


tools_mod = Extension('MOODS._tools',
                           sources=['../core/tools_wrap.cxx', '../core/moods_tools.cpp'],
                           include_dirs=["../core/"],
                           extra_compile_args=['-march=native', '-O3', '-fPIC', '--std=c++0x'],
                           )

setup (name = 'MOODS',
       version = '1.9',
       description = 'MOODS: Motif Occurrence Detection Suite',
       ext_modules = [tools_mod],
       py_modules = ["MOODS.tools"],
)
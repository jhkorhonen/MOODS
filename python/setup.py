from distutils.core import setup, Extension


module1 = Extension('MOODS._cmodule',
                    sources = ['cport.cc'],
                    include_dirs=["../src"],
                    library_dirs=["../src"],
                    libraries=['pssm'])

setup (name = 'MOODS',
       version = '1.0.2',
       description = 'MOODS: Motif Occurrence Detection Suite',
       packages = ['MOODS',],
       ext_modules = [module1])
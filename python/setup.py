"""
setup.py file for MOODS
"""

from setuptools import setup, Extension
from pathlib import Path
import subprocess
import shutil
import sys

here = Path(__file__).parent.resolve()

# Locate C++ code under core/ ; this folder should exist within the python/ folder or one level above it
core_dir = here / "core"
if not core_dir.exists():
    core_dir = (here.parent / "core")
if not core_dir.exists():
    raise RuntimeError(f"core/ not found next to or above {here}")

moods_pkg_dir = here / "MOODS"
moods_pkg_dir.mkdir(exist_ok=True)

def run_swig():
    """Generate Python wrappers and C++ bindings with swig before setup.py scans directory."""
    if not shutil.which("swig"):
        raise RuntimeError("swig not found in PATH; install swig and retry")

    # Generate wrappers for these interface files
    for mod in ["scan", "tools", "misc", "parsers"]:
        i_file = core_dir / f"{mod}.i"
        py_out = moods_pkg_dir / f"{mod}.py"
        cxx_out = core_dir / f"{mod}_wrap.cxx"

        # Regenerate if missing
        # Note: You can force regen by deleting either output file.
        if not py_out.exists() or not cxx_out.exists():
            cmd = [
                "swig",
                "-c++",
                "-python",
                "-outdir", str(moods_pkg_dir),
                str(i_file),
            ]
            subprocess.check_call(cmd, cwd=here)

def compile_args():
    if sys.platform.startswith("win"):
        # MSVC: keep minimal, portable flags.
        return ["/O2", "/EHsc"]
    else:
        # GCC/Clang
        return ["-O3", "-std=c++11", "-fPIC"]

# Ensure wrappers exist before setuptools scans packages/files
run_swig()

common_includes = [str(core_dir)]
common_compile_args = compile_args()

tools_mod = Extension(
    "MOODS._tools",
    sources=[
        str(core_dir / "tools_wrap.cxx"),
        str(core_dir / "moods_tools.cpp"),
        str(core_dir / "moods_misc.cpp"),
        str(core_dir / "match_types.cpp"),
    ],
    include_dirs=common_includes,
    extra_compile_args=common_compile_args,
)

scan_mod = Extension(
    "MOODS._scan",
    sources=[
        str(core_dir / "scan_wrap.cxx"),
        str(core_dir / "moods_scan.cpp"),
        str(core_dir / "motif_0.cpp"),
        str(core_dir / "motif_h.cpp"),
        str(core_dir / "moods_misc.cpp"),
        str(core_dir / "scanner.cpp"),
        str(core_dir / "moods_tools.cpp"),
        str(core_dir / "match_types.cpp"),
    ],
    include_dirs=common_includes,
    extra_compile_args=common_compile_args,
)

parsers_mod = Extension(
    "MOODS._parsers",
    sources=[
        str(core_dir / "parsers_wrap.cxx"),
        str(core_dir / "moods_parsers.cpp"),
        str(core_dir / "moods_misc.cpp"),
        str(core_dir / "moods_tools.cpp"),
        str(core_dir / "match_types.cpp"),
    ],
    include_dirs=common_includes,
    extra_compile_args=common_compile_args,
)

readme_file = here / "readme.MD"
if readme_file.exists():
    with open(here / 'readme.MD') as f:
        long_description = f.read()
else:
    long_description = "MOODS: Motif Occurrence Detection Suite (Python bindings)"

setup(
    name="MOODS-python",
    version="1.9.4.1",
    description="MOODS: Motif Occurrence Detection Suite",
    long_description=long_description,
    long_description_content_type="text/markdown",
    maintainer="Janne H. Korhonen",
    maintainer_email="janne.h.korhonen@gmail.com",
    url="https://www.cs.helsinki.fi/group/pssmfind/",
    license="GPLv3 / Biopython license",
    ext_modules=[tools_mod, scan_mod, parsers_mod],    
    packages=["MOODS"],
    scripts=['scripts/moods-dna.py'],
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords="PWM, PSSM, motif scan",
)

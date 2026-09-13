"""Build the schubmult_cpp extension (Cython + the C++ kernels in this directory).

    cd cpp && python setup_ext.py build_ext --inplace

Requires Cython, the SymEngine C++ headers/library (defaults to the active conda env,
override with SYMENGINE_PREFIX), and an importable symengine.py checkout with its .pxd
files (as installed in the schubmult_312 env). The resulting schubmult_cpp*.so is placed
in this directory; put it on sys.path (or copy it into src/schubmult/) to use it.
"""

import os
import sys
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup

import symengine

here = Path(__file__).parent.resolve()
prefix = Path(os.environ.get("SYMENGINE_PREFIX", os.environ.get("CONDA_PREFIX", sys.prefix)))
maxn = os.environ.get("MAXN", "32")

# symengine.py ships its Cython declarations next to the wrapper module
se_pkg = Path(symengine.__file__).parent
se_lib = se_pkg / "lib"
if not (se_lib / "symengine_wrapper.pxd").exists():
    sys.exit(f"symengine_wrapper.pxd not found in {se_lib}; need a symengine.py source install")

ext = Extension(
    "schubmult_cpp",
    sources=[str(here / "schubmult_cpp.pyx")],
    language="c++",
    include_dirs=[str(here), str(prefix / "include")],
    library_dirs=[str(prefix / "lib")],
    libraries=["symengine", "flint", "gmp"],
    runtime_library_dirs=[str(prefix / "lib")],
    extra_compile_args=["-O3", "-march=native", "-std=c++17", f"-DMAXN={maxn}"],
)

setup(
    name="schubmult_cpp",
    ext_modules=cythonize([ext], include_path=[str(se_lib), str(se_pkg.parent)], compiler_directives={"language_level": 3}),
    script_args=sys.argv[1:] or ["build_ext", "--inplace"],
)

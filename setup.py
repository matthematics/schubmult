"""Builds the schubmult_cpp extension: the C++ multiplication kernels in cpp/, exposed through the
plain Python C API (cpp/schubmult_module.cpp). Package metadata lives in pyproject.toml.

The extension needs only a C++17 compiler and the Python headers. Coefficients are handled as
Python `symengine` objects, so no SymEngine C++ headers or libraries are involved and the
runtime `symengine` wheel from PyPI is all that is required.

Environment knobs: MAXN (default 32, largest permutation size), SCHUBMULT_NATIVE=1 (tune for the
build machine's CPU; off by default so wheels are portable).
"""

import os
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

here = Path(__file__).parent.resolve()
maxn = os.environ.get("MAXN", "32")


class BuildExt(build_ext):
    """Per-compiler flags (MSVC spells them differently)."""

    def build_extensions(self):
        msvc = self.compiler.compiler_type == "msvc"
        for ext in self.extensions:
            if msvc:
                ext.extra_compile_args = ["/std:c++17", "/O2", "/EHsc"]
            else:
                ext.extra_compile_args = ["-std=c++17", "-O3"]
                if os.environ.get("SCHUBMULT_NATIVE"):
                    ext.extra_compile_args.append("-march=native")
                if sys.platform == "darwin":
                    ext.extra_compile_args.append("-mmacosx-version-min=10.14")  # std::filesystem-free, but needs a C++17 libc++
        super().build_extensions()


ext = Extension(
    "schubmult.schubmult_cpp",
    sources=[str(here / "cpp" / "schubmult_module.cpp")],
    depends=sorted(str(p) for p in (here / "cpp").glob("*.h")),  # setuptools only tracks sources otherwise
    include_dirs=[str(here / "cpp")],
    define_macros=[("MAXN", maxn), ("SCHUB_PYEXPR", "1")],
    language="c++",
)

setup(ext_modules=[ext], cmdclass={"build_ext": BuildExt})

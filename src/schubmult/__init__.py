from __future__ import annotations

import ast
import importlib
import inspect
import pkgutil
import warnings
from pathlib import Path
from typing import Dict, List, Set

__version__ = "4.0.0dev"

"""
Top-level schubmult package initializer.

Dynamically import and re-export public classes defined in submodules
under this package. Imports are best-effort: modules that fail to import
are skipped so that a top-level `import schubmult` does not raise due to
optional heavy dependencies (e.g. pulp/mosek).
"""

__all__: List[str] = []
# map exported name -> module that originally defined it (for lazy re-importing)
_module_map: Dict[str, str] = {}

# Build a mapping name -> module (relative import path) by scanning source files.
# This uses AST only (no module imports), so it's safe at package import time.
_package_root = Path(__file__).resolve().parent
_skipped_dirs = {"__pycache__", "tests", "docs", "scripts", "build"}

for py in _package_root.rglob("*.py"):
    # skip this file and any in skipped directories
    if py.samefile(Path(__file__)):
        continue
    if any(part in _skipped_dirs for part in py.parts):
        continue
    # compute module name like "schubmult.sub.module"
    try:
        rel = py.relative_to(_package_root.parent)
    except Exception:
        continue
    modname = ".".join(rel.with_suffix("").parts)
    # parse AST and collect top-level defs/assigns (public names)
    try:
        src = py.read_text(encoding="utf8")
        tree = ast.parse(src)
    except Exception:
        continue
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = node.name
            if not name.startswith("_") and name not in _module_map:
                _module_map[name] = modname
        elif isinstance(node, ast.Assign):
            # collect simple targets (NAME = ...)
            for target in node.targets:
                if isinstance(target, ast.Name):
                    name = target.id
                    if not name.startswith("_") and name not in _module_map:
                        _module_map[name] = modname

# Export the discovered public names; keep as a sorted tuple for reproducibility.
__all__ = tuple(sorted(set(_module_map.keys())))

# Commonly exported helpers that we prefer to lazy-load
# we do not import them eagerly to avoid triggering heavy dependencies.
_lazy_exports = {
    "Sx": "schubmult.rings.schubert_ring",
    "uncode": "schubmult.schub_lib.perm_lib",
}
for name in _lazy_exports:
    if name not in __all__:
        __all__.append(name)

# finalize __all__
__all__ = tuple(sorted(set(__all__)))

def __getattr__(name: str):
    """
    Lazy import exported names that were not successfully imported at package
    initialization time. This allows `import schubmult` to succeed even if
    some optional dependencies (used by particular submodules) are missing.
    """
    if name in globals():
        return globals()[name]
    # known lazy exports
    if name in _lazy_exports:
        try:
            mod = importlib.import_module(_lazy_exports[name])
            val = getattr(mod, name)
            globals()[name] = val
            return val
        except Exception as e:
            raise AttributeError(f"cannot import {name!r} from {_lazy_exports[name]!r}: {e}") from e
    # names discovered during the initial scan
    if name in _module_map:
        try:
            mod = importlib.import_module(_module_map[name])
            val = getattr(mod, name)
            globals()[name] = val
            return val
        except Exception as e:
            raise AttributeError(f"cannot import {name!r} from {_module_map[name]!r}: {e}") from e
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    # include lazily-exported names in dir()
    return sorted(set(list(globals().keys()) + list(__all__)))

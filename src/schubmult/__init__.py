from __future__ import annotations

import ast
import importlib
from pathlib import Path
from typing import Dict, List

# Version is managed by setuptools-scm
try:
    from schubmult._version import version as __version__
except ImportError:
    __version__ = "4.0.0.dev0"

"""
Top-level schubmult package initializer.

Lazily import and re-export public classes defined in submodules.
Uses deferred AST scanning - only scans files when __dir__() or an unknown
attribute is accessed for the first time.
"""

# Commonly exported helpers - only these are in __all__ initially
_lazy_exports = {
    # Core Schubert structures
    "Permutation": "schubmult.schub_lib.permutation",
    "uncode": "schubmult.schub_lib.permutation",
    "RCGraph": "schubmult.schub_lib.rc_graph",
    "ChuteMoveElement": "schubmult.schub_lib.chute_move_element",
    "ChuteMoveRing": "schubmult.rings.combinatorial.chute_move_ring",
    "BPD": "schubmult.schub_lib.bpd",
    "TileType": "schubmult.schub_lib.bpd",
    "HPD": "schubmult.schub_lib.hpd",
    "HPDTile": "schubmult.schub_lib.hpd",
    # Rings
    "Sx": "schubmult.rings.schubert.schubert_ring",
    "DSx": "schubmult.rings.schubert.schubert_ring",
    "ASx": "schubmult.rings.free_algebra",
    "ADSx": "schubmult.rings.free_algebra",
    "FA": "schubmult.rings.free_algebra",
    "DoubleSchubertRing": "schubmult.rings.schubert.schubert_ring",
    "SingleSchubertRing": "schubmult.rings.schubert.schubert_ring",
    "DoubleSchubertElement": "schubmult.rings.schubert.schubert_ring",
    "QSym": "schubmult.rings.quasisymmetric_functions",
    "QSymElement": "schubmult.rings.quasisymmetric_functions",
    "RCGraphRing": "schubmult.rings.combinatorial.rc_graph_ring",
    "GrassRCGraphRing": "schubmult.rings.combinatorial.rc_graph_ring",
    "RCGraphRingElement": "schubmult.rings.combinatorial.rc_graph_ring",
    "EGPlacticRing": "schubmult.rings.combinatorial.eg_plactic_ring",
    "EGPlacticRingElement": "schubmult.rings.combinatorial.eg_plactic_ring",
    "BPDRing": "schubmult.rings.combinatorial.bpd_ring",
    "BPDRingElement": "schubmult.rings.combinatorial.bpd_ring",
    "SchubertMonomialRing": "schubmult.rings.combinatorial.schubert_monomial_ring",
    "SchubertMonomialRingElement": "schubmult.rings.combinatorial.schubert_monomial_ring",
    "PlacticAlgebra": "schubmult.rings.combinatorial.plactic_algebra",
    "PlacticAlgebraElement": "schubmult.rings.combinatorial.plactic_algebra",
    "CoxeterKnuthRing": "schubmult.rings.ck_ring",
    "FreeAlgebra": "schubmult.rings.free_algebra",
    "NSym": "schubmult.rings.nsym",
    "NilHeckeRing": "schubmult.rings.nil_hecke",
    "PolynomialAlgebra": "schubmult.rings.polynomial_algebra",
    "TensorRing": "schubmult.rings.tensor_ring",
    "BaseRing": "schubmult.rings.base_ring",
    "BaseRingElement": "schubmult.rings.base_ring",
    # Variables
    "GeneratingSet": "schubmult.symbolic.poly.variables",
    "pull_out_var": "schubmult.utils.schub_lib",
    # Crystal structures
    "RootTableau": "schubmult.schub_lib.root_tableau",
    "Plactic": "schubmult.schub_lib.plactic",
    "NilPlactic": "schubmult.schub_lib.nilplactic",
    "NilPlacticAlgebra": "schubmult.rings.combinatorial.plactic_algebra",
    # Free algebra bases
    "FreeAlgebraBasis": "schubmult.rings.free_algebra",
    "WordBasis": "schubmult.rings.free_algebra",
    "JBasis": "schubmult.rings.free_algebra",
    "JTBasis": "schubmult.rings.free_algebra",
    "SchubertBasis": "schubmult.rings.free_algebra",
    "SchubertSchurBasis": "schubmult.rings.free_algebra",
    "ElementaryBasis": "schubmult.rings.free_algebra",
    "NElementaryBasis": "schubmult.rings.free_algebra",
    "ZBasis": "schubmult.rings.free_algebra",
}

__all__: List[str] = list(_lazy_exports.keys())

# Cache for lazy discovery - populated on first access
_module_map: Dict[str, str] = {}
_scanned = False
_package_root = Path(__file__).resolve().parent

_package_root = Path(__file__).resolve().parent

def _scan_modules():
    """Deferred AST scan - only called when needed"""
    global _scanned  # noqa: PLW0603
    if _scanned:
        return
    _scanned = True

    _skipped_dirs = {"__pycache__", "tests", "docs", "scripts", "build", "_scripts"}

    for py in _package_root.rglob("*.py"):
        if py.samefile(Path(__file__)):
            continue
        if any(part in _skipped_dirs for part in py.parts):
            continue
        try:
            rel = py.relative_to(_package_root.parent)
        except Exception:
            continue
        modname = ".".join(rel.with_suffix("").parts)

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
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        name = target.id
                        if not name.startswith("_") and name not in _module_map:
                            _module_map[name] = modname

def __getattr__(name: str):
    """
    Lazy import exported names. This allows `import schubmult` to succeed
    even if some optional dependencies are missing.
    """
    if name in globals():
        return globals()[name]

    # Known lazy exports - try these first without scanning
    if name in _lazy_exports:
        try:
            mod = importlib.import_module(_lazy_exports[name])
            val = getattr(mod, name)
            globals()[name] = val
            return val
        except Exception as e:
            raise AttributeError(f"cannot import {name!r} from {_lazy_exports[name]!r}: {e}") from e

    # Not in known exports - scan modules if we haven't yet
    _scan_modules()

    if name in _module_map:
        modname = _module_map[name]
        try:
            mod = importlib.import_module(modname)
            val = getattr(mod, name)
            globals()[name] = val
            return val
        except Exception as e:
            # If the failing module lives under schubmult._scripts, be lenient:
            # warn the user but return a lightweight proxy that raises a clear
            # error when used. This lets `import schubmult` succeed while still
            # giving a helpful message if the script is actually executed.
            if modname.startswith(__name__ + "._scripts"):
                def _broken_callable(*_, _name=name, _mod=modname, _err=e, **__):
                    raise RuntimeError(
                        f"attempted to use {_name!r} from {_mod!r} but importing the module failed: {_err!r}",
                    )

                # attach some attributes to help introspection / debugging
                broken = _broken_callable
                broken.__name__ = f"broken_{name}"
                broken.__doc__ = f"Placeholder for {name} from {modname}; import failed: {e!r}"
                globals()[name] = broken
                return broken

            # otherwise behave as before
            raise AttributeError(f"cannot import {name!r} from {_module_map[name]!r}: {e}") from e

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    """Include lazily-exported names in dir()"""
    _scan_modules()  # Ensure we have full module map
    return sorted(set(list(globals().keys()) + list(_lazy_exports.keys()) + list(_module_map.keys())))

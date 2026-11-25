from __future__ import annotations

import ast
import importlib
import importlib.util
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List

if TYPE_CHECKING:
    from schubmult.combinatorics.anti_rc_graph import (
        AntiRCGraph,
    )
    from schubmult.combinatorics.bpd import (
        BPD,
        TileType,
    )
    from schubmult.combinatorics.crystal_graph import (
        CrystalGraphTensor,
    )
    from schubmult.combinatorics.increasing_tableau import (
        IncreasingTableau,
    )
    from schubmult.combinatorics.nilplactic import (
        NilPlactic,
    )
    from schubmult.combinatorics.permutation import (
        Permutation,
        uncode,
    )
    from schubmult.combinatorics.pipe_dream import (
        PipeDream,
    )
    from schubmult.combinatorics.plactic import (
        Plactic,
    )
    from schubmult.combinatorics.rc_graph import (
        RCGraph,
    )
    from schubmult.combinatorics.set_valued_tableau import (
        SetValuedTableau,
    )
    from schubmult.combinatorics.wc_graph import (
        WCGraph,
    )
    from schubmult.rings.base_ring import (
        BaseRing,
        BaseRingElement,
    )
    from schubmult.rings.combinatorial.bounded_rc_factor_algebra import (
        BoundedRCFactorAlgebra,
    )
    from schubmult.rings.combinatorial.bounded_wc_factor_algebra import (
        BoundedWCFactorAlgebra,
    )
    from schubmult.rings.combinatorial.bpd_ring import (
        BPDRing,
        BPDRingElement,
    )
    from schubmult.rings.combinatorial.chute_move_ring import (
        ChuteMoveRing,
    )
    from schubmult.rings.combinatorial.dual_rc_graph_ring import (
        DualRCGraphRing,
    )
    from schubmult.rings.combinatorial.eg_plactic_ring import (
        EGPlacticRing,
        EGPlacticRingElement,
    )
    from schubmult.rings.combinatorial.forest_rc_ring import (
        ForestRCGraphRing,
    )
    from schubmult.rings.combinatorial.grass_tensor_algebra import (
        GrassTensorAlgebra,
        GrassTensorAlgebraElement,
    )
    from schubmult.rings.combinatorial.plactic_algebra import (
        NilPlacticAlgebra,
        PlacticAlgebra,
        PlacticAlgebraElement,
    )
    from schubmult.rings.combinatorial.rc_graph_ring import (
        GrassRCGraphRing,
        RCGraphRing,
        RCGraphRingElement,
    )
    from schubmult.rings.combinatorial.schubert_monomial_ring import (
        SchubertMonomialRing,
        SchubertMonomialRingElement,
    )
    from schubmult.rings.combinatorial.wc_graph_ring import (
        WCGraphRing,
        WCGraphRingElement,
    )
    from schubmult.rings.direct_product_ring import (
        DirectProductRing,
        DirectProductRingElement,
    )
    from schubmult.rings.free_algebra import (
        FA,
        ADSx,
        AGx,
        ASx,
        CompositionSchubertBasis,
        ElementaryBasis,
        FreeAlgebra,
        FreeAlgebraBasis,
        GrothendieckBasis,
        JBasis,
        JTBasis,
        NElementaryBasis,
        SchubertBasis,
        SchubertSchurBasis,
        WordBasis,
        ZBasis,
    )
    from schubmult.rings.free_algebra.forest_basis import (
        ForestBasis,
        ForestDual,
    )
    from schubmult.rings.nsym import (
        NSym,
    )
    from schubmult.rings.polynomial_algebra import (
        PolynomialAlgebra,
    )
    from schubmult.rings.quasisymmetric_functions import (
        QSym,
        QSymElement,
    )
    from schubmult.rings.schubert.double_grothendieck_ring import (
        DGx,
        DoubleGrothendieckElement,
        DoubleGrothendieckRing,
    )
    from schubmult.rings.schubert.grothendieck_ring import (
        GrothendieckElement,
        GrothendieckRing,
        Gx,
    )
    from schubmult.rings.schubert.nil_hecke import (
        NilHeckeRing,
    )
    from schubmult.rings.schubert.quantum_double_grothendieck_ring import (
        QDGx,
        QuantumDoubleGrothendieckElement,
        QuantumDoubleGrothendieckRing,
    )
    from schubmult.rings.schubert.schubert_ring import (
        DoubleSchubertElement,
        DoubleSchubertRing,
        DSx,
        SingleSchubertRing,
        Sx,
    )
    from schubmult.rings.tensor_ring import (
        TensorRing,
    )
    from schubmult.rings.thompson_algebra import (
        ThompsonAlgebra,
    )
    from schubmult.symbolic.poly.variables import (
        ChangeOfVariables,
        GeneratingSet,
    )
    from schubmult.utils.schub_lib import (
        pull_out_var,
    )

# Version is managed by setuptools-scm
try:
    from schubmult._version import version as __version__
except ImportError:
    __version__ = "5.1.0"

"""
Top-level schubmult package initializer.

Lazily import and re-export public classes defined in non-experimental submodules.
Uses deferred AST scanning - only scans files when __dir__() or an unknown
attribute is accessed for the first time.
"""

# Commonly exported helpers - only these are in __all__ initially
_lazy_exports = {
    # Core Schubert structures
    "Permutation": "schubmult.combinatorics.permutation",
    "uncode": "schubmult.combinatorics.permutation",
    "RCGraph": "schubmult.combinatorics.rc_graph",
    "WCGraph": "schubmult.combinatorics.wc_graph",
    "AntiRCGraph": "schubmult.combinatorics.anti_rc_graph",
    "BoundedRCFactorAlgebra": "schubmult.rings.combinatorial.bounded_rc_factor_algebra",
    "BoundedWCFactorAlgebra": "schubmult.rings.combinatorial.bounded_wc_factor_algebra",
    "ThompsonAlgebra": "schubmult.rings.thompson_algebra",
    "ChuteMoveRing": "schubmult.rings.combinatorial.chute_move_ring",
    "BPD": "schubmult.combinatorics.bpd",
    "ForestRCGraphRing": "schubmult.rings.combinatorial.forest_rc_ring",
    "TileType": "schubmult.combinatorics.bpd",
    # Rings
    "Gx": "schubmult.rings.schubert.grothendieck_ring",
    "DGx": "schubmult.rings.schubert.double_grothendieck_ring",
    "DoubleGrothendieckRing": "schubmult.rings.schubert.double_grothendieck_ring",
    "DoubleGrothendieckElement": "schubmult.rings.schubert.double_grothendieck_ring",
    "QDGx": "schubmult.rings.schubert.quantum_double_grothendieck_ring",
    "QuantumDoubleGrothendieckRing": "schubmult.rings.schubert.quantum_double_grothendieck_ring",
    "QuantumDoubleGrothendieckElement": "schubmult.rings.schubert.quantum_double_grothendieck_ring",
    "GrothendieckRing": "schubmult.rings.schubert.grothendieck_ring",
    "GrothendieckElement": "schubmult.rings.schubert.grothendieck_ring",
    "Sx": "schubmult.rings.schubert.schubert_ring",
    "DSx": "schubmult.rings.schubert.schubert_ring",
    "ASx": "schubmult.rings.free_algebra",
    "AGx": "schubmult.rings.free_algebra",
    "ADSx": "schubmult.rings.free_algebra",
    "FA": "schubmult.rings.free_algebra",
    "PipeDream": "schubmult.combinatorics.pipe_dream",
    "IncreasingTableau": "schubmult.combinatorics.increasing_tableau",
    "SetValuedTableau": "schubmult.combinatorics.set_valued_tableau",
    "DualRCGraphRing": "schubmult.rings.combinatorial.dual_rc_graph_ring",
    "ForestDual": "schubmult.rings.free_algebra.forest_basis",
    "ForestBasis": "schubmult.rings.free_algebra.forest_basis",
    "DoubleSchubertRing": "schubmult.rings.schubert.schubert_ring",
    "SingleSchubertRing": "schubmult.rings.schubert.schubert_ring",
    "DoubleSchubertElement": "schubmult.rings.schubert.schubert_ring",
    "QSym": "schubmult.rings.quasisymmetric_functions",
    "QSymElement": "schubmult.rings.quasisymmetric_functions",
    "RCGraphRing": "schubmult.rings.combinatorial.rc_graph_ring",
    "GrassRCGraphRing": "schubmult.rings.combinatorial.rc_graph_ring",
    "GrassTensorAlgebra": "schubmult.rings.combinatorial.grass_tensor_algebra",
    "GrassTensorAlgebraElement": "schubmult.rings.combinatorial.grass_tensor_algebra",
    "RCGraphRingElement": "schubmult.rings.combinatorial.rc_graph_ring",
    "WCGraphRing": "schubmult.rings.combinatorial.wc_graph_ring",
    "WCGraphRingElement": "schubmult.rings.combinatorial.wc_graph_ring",
    "EGPlacticRing": "schubmult.rings.combinatorial.eg_plactic_ring",
    "EGPlacticRingElement": "schubmult.rings.combinatorial.eg_plactic_ring",
    "BPDRing": "schubmult.rings.combinatorial.bpd_ring",
    "BPDRingElement": "schubmult.rings.combinatorial.bpd_ring",
    "CrystalGraphTensor": "schubmult.combinatorics.crystal_graph",
    "SchubertMonomialRing": "schubmult.rings.combinatorial.schubert_monomial_ring",
    "SchubertMonomialRingElement": "schubmult.rings.combinatorial.schubert_monomial_ring",
    "PlacticAlgebra": "schubmult.rings.combinatorial.plactic_algebra",
    "PlacticAlgebraElement": "schubmult.rings.combinatorial.plactic_algebra",
    "FreeAlgebra": "schubmult.rings.free_algebra",
    "NSym": "schubmult.rings.nsym",
    "NilHeckeRing": "schubmult.rings.schubert.nil_hecke",
    "PolynomialAlgebra": "schubmult.rings.polynomial_algebra",
    "TensorRing": "schubmult.rings.tensor_ring",
    "DirectProductRing": "schubmult.rings.direct_product_ring",
    "DirectProductRingElement": "schubmult.rings.direct_product_ring",
    "BaseRing": "schubmult.rings.base_ring",
    "BaseRingElement": "schubmult.rings.base_ring",
    # Variables
    "GeneratingSet": "schubmult.symbolic.poly.variables",
    "ChangeOfVariables": "schubmult.symbolic.poly.variables",
    "pull_out_var": "schubmult.utils.schub_lib",
    # Crystal structures
    "Plactic": "schubmult.combinatorics.plactic",
    "NilPlactic": "schubmult.combinatorics.nilplactic",
    "NilPlacticAlgebra": "schubmult.rings.combinatorial.plactic_algebra",
    # Free algebra bases
    "FreeAlgebraBasis": "schubmult.rings.free_algebra",
    "WordBasis": "schubmult.rings.free_algebra",
    "JBasis": "schubmult.rings.free_algebra",
    "JTBasis": "schubmult.rings.free_algebra",
    "SchubertBasis": "schubmult.rings.free_algebra",
    "CompositionSchubertBasis": "schubmult.rings.free_algebra",
    "SchubertSchurBasis": "schubmult.rings.free_algebra",
    "ElementaryBasis": "schubmult.rings.free_algebra",
    "NElementaryBasis": "schubmult.rings.free_algebra",
    "GrothendieckBasis": "schubmult.rings.free_algebra",
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
    experimental_root = _package_root / "combinatorics" / "experimental"

    for py in _package_root.rglob("*.py"):
        if py.samefile(Path(__file__)):
            continue
        if experimental_root in py.parents:
            continue
        if any(part in _skipped_dirs for part in py.parts):
            continue
        if py.stem.startswith("_") and py.stem != "__init__":
            continue  # private modules (e.g. mult._accel) may deliberately reuse public names
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
            if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                name = node.name
                if not name.startswith("_") and name not in _module_map:
                    _module_map[name] = modname
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        name = target.id
                        if not name.startswith("_") and name not in _module_map:
                            _module_map[name] = modname

def __getattr__(name: str) -> Any:
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

    # Submodules (e.g. ``from schubmult import schubmult_cpp``) must not trigger the slow scan
    if importlib.util.find_spec(f"{__name__}.{name}") is not None:
        return importlib.import_module(f"{__name__}.{name}")

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

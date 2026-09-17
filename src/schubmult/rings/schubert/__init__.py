"""Schubert-family rings: the main user-facing algebra interface.

The most commonly used entry points are the ring *instances*:

- ``Sx``: ordinary (single) Schubert polynomials, ``Sx([3, 1, 2]) * Sx([2, 1, 3])``.
- ``DSx``: double Schubert polynomials (second alphabet ``y``).
- ``Gx`` / ``DGx``: (double) Grothendieck polynomials.
- ``QSx`` / ``QDSx``: quantum (double) Schubert polynomials.
- ``QPSx`` / ``QPDSx``: parabolic quantum (double) Schubert polynomials.

Each instance is an object of the corresponding ``*Ring`` class; calling it with
a permutation (or Lehmer code, or a polynomial expression) yields a ``*Element``
that supports ``+``, ``*``, ``.expand()``, and conversion between bases. All ring
classes derive from `BaseSchubertRing` and dispatch their products to the kernels
in `schubmult.mult`.

Everything here is imported lazily (see ``__getattr__``) so ``import schubmult``
stays fast.
"""

from __future__ import annotations

import importlib

_lazy_exports = {
    "BaseSchubertElement": "schubmult.rings.schubert.base_schubert_ring",
    "BaseSchubertRing": "schubmult.rings.schubert.base_schubert_ring",
    "DGx": "schubmult.rings.schubert.double_grothendieck_ring",
    "DSx": "schubmult.rings.schubert.double_schubert_ring",
    "Gx": "schubmult.rings.schubert.grothendieck_ring",
    "GrothendieckElement": "schubmult.rings.schubert.grothendieck_ring",
    "GrothendieckRing": "schubmult.rings.schubert.grothendieck_ring",
    "DoubleGrothendieckElement": "schubmult.rings.schubert.double_grothendieck_ring",
    "DoubleGrothendieckRing": "schubmult.rings.schubert.double_grothendieck_ring",
    "DoubleSchubertElement": "schubmult.rings.schubert.double_schubert_ring",
    "DoubleSchubertRing": "schubmult.rings.schubert.double_schubert_ring",
    "ElemDoubleSchubertRing": "schubmult.rings.schubert.double_schubert_ring",
    "NilHeckeElement": "schubmult.rings.schubert.nil_hecke",
    "NilHeckeRing": "schubmult.rings.schubert.nil_hecke",
    "ParabolicQuantumDoubleSchubertElement": "schubmult.rings.schubert.parabolic_quantum_double_schubert_ring",
    "ParabolicQuantumDoubleSchubertRing": "schubmult.rings.schubert.parabolic_quantum_double_schubert_ring",
    "QDSx": "schubmult.rings.schubert.quantum_double_schubert_ring",
    "QPDSx": "schubmult.rings.schubert.parabolic_quantum_double_schubert_ring",
    "QPSx": "schubmult.rings.schubert.parabolic_quantum_schubert_ring",
    "QSx": "schubmult.rings.schubert.quantum_schubert_ring",
    "QuantumDoubleSchubertElement": "schubmult.rings.schubert.quantum_double_schubert_ring",
    "QuantumDoubleSchubertRing": "schubmult.rings.schubert.quantum_double_schubert_ring",
    "QuantumSingleSchubertRing": "schubmult.rings.schubert.quantum_schubert_ring",
    "SeparatedDescentsRing": "schubmult.rings.schubert.separated_descents",
    "SeparatedDescentsRingElement": "schubmult.rings.schubert.separated_descents",
    "SingleSchubertRing": "schubmult.rings.schubert.schubert_ring",
    "Sx": "schubmult.rings.schubert.schubert_ring",
    "make_parabolic_quantum_basis": "schubmult.rings.schubert.parabolic_quantum_double_schubert_ring",
}

__all__ = [
    "BaseSchubertElement",
    "BaseSchubertRing",
    "DGx",
    "DSx",
    "DoubleGrothendieckElement",
    "DoubleGrothendieckRing",
    "DoubleSchubertElement",
    "DoubleSchubertRing",
    "ElemDoubleSchubertRing",
    "GrothendieckElement",
    "GrothendieckRing",
    "Gx",
    "NilHeckeElement",
    "NilHeckeRing",
    "ParabolicQuantumDoubleSchubertElement",
    "ParabolicQuantumDoubleSchubertRing",
    "QDSx",
    "QPDSx",
    "QPSx",
    "QSx",
    "QuantumDoubleSchubertElement",
    "QuantumDoubleSchubertRing",
    "QuantumSingleSchubertRing",
    "SeparatedDescentsRing",
    "SeparatedDescentsRingElement",
    "SingleSchubertRing",
    "Sx",
    "make_parabolic_quantum_basis",
]


def __getattr__(name: str):
    """Lazily import and cache the requested export from its defining submodule."""
    if name not in _lazy_exports:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    mod = importlib.import_module(_lazy_exports[name])
    val = getattr(mod, name)
    globals()[name] = val
    return val


def __dir__():
    """Include lazily-exported names in ``dir()``."""
    return sorted(set(list(globals().keys()) + list(__all__)))

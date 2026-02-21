from __future__ import annotations

import importlib

_lazy_exports = {
    "BaseSchubertElement": "schubmult.rings.schubert.base_schubert_ring",
    "BaseSchubertRing": "schubmult.rings.schubert.base_schubert_ring",
    "DSx": "schubmult.rings.schubert.double_schubert_ring",
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
    "DSx",
    "DoubleSchubertElement",
    "DoubleSchubertRing",
    "ElemDoubleSchubertRing",
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
    if name not in _lazy_exports:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    mod = importlib.import_module(_lazy_exports[name])
    val = getattr(mod, name)
    globals()[name] = val
    return val


def __dir__():
    return sorted(set(list(globals().keys()) + list(__all__)))

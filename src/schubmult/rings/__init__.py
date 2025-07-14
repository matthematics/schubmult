from .free_algebra import FA, FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import FreeAlgebraBasis, SchubertBasis, SchubertSchurBasis, SeparatedDescentsBasis, WordBasis
from .nil_hecke import NilHeckeRing
from .quantum_schubert_ring import QDSx, QPDSx, QPSx, QSx, QuantumDoubleSchubertElement, QuantumDoubleSchubertRing, QuantumSingleSchubertRing, make_parabolic_quantum_basis
from .schubert_ring import DoubleSchubertElement, DoubleSchubertRing, DSx, ElemDoubleSchubertRing, SingleSchubertRing, Sx
from .separated_descents import SeparatedDescentsRing, SeparatedDescentsRingElement
from .tensor_ring import TensorRing, TensorRingElement

# ASx = #SeparatedDescentsRing(Sx([]).ring)
ASx = FreeAlgebra(basis=SchubertBasis)
ADSx = SeparatedDescentsRing(DSx([]).ring)

__all__ = [
    "FA",
    "ADSx",
    "ASx",
    "DSx",
    "DoubleSchubertElement",
    "DoubleSchubertRing",
    "ElemDoubleSchubertRing",
    "FreeAlgebra",
    "FreeAlgebraBasis",
    "FreeAlgebraElement",
    "NilHeckeRing",
    "QDSx",
    "QPDSx",
    "QPSx",
    "QSx",
    "QuantumDoubleSchubertElement",
    "QuantumDoubleSchubertRing",
    "QuantumSingleSchubertRing",
    "SchubertSchurBasis",
    "SeparatedDescentsBasis",
    "SeparatedDescentsRing",
    "SeparatedDescentsRingElement",
    "SingleSchubertRing",
    "Sx",
    "TensorRing",
    "TensorRingElement",
    "WordBasis",
    "make_parabolic_quantum_basis",
]

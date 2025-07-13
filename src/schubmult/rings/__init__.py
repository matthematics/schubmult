from .free_algebra import FA, FreeAlgebra, FreeAlgebraElement, SchubertBasis, WordBasis
from .nil_hecke import NilHeckeRing
from .quantum_schubert_ring import QDSx, QPDSx, QPSx, QSx, QuantumDoubleSchubertElement, QuantumDoubleSchubertRing, QuantumSingleSchubertRing, make_parabolic_quantum_basis
from .schubert_ring import DoubleSchubertElement, DoubleSchubertRing, DSx, ElemDoubleSchubertRing, SingleSchubertRing, Sx
from .separated_descents import SeparatedDescentsRing, SeparatedDescentsRingElement
from .tensor_ring import TensorRing, TensorRingElement

#ASx = #SeparatedDescentsRing(Sx([]).ring)
ASx = FreeAlgebra(basis=SchubertBasis())
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
    "FreeAlgebraElement",
    "NilHeckeRing",
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
    "TensorRing",
    "TensorRingElement",
    "make_parabolic_quantum_basis",
]

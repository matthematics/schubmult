from .free_algebra import FA, FreeAlgebra, FreeAlgebraElement
from .free_algebra_basis import ElementaryBasis, FreeAlgebraBasis, JBasis, NElementaryBasis, SchubertBasis, SchubertSchurBasis, SeparatedDescentsBasis, WordBasis, ZBasis
from .nil_hecke import NilHeckeRing
from .polynomial_algebra import PA, PolynomialAlgebra, PolynomialAlgebraElement
from .polynomial_basis import ElemSymPolyBasis, MonomialBasis, PolynomialBasis, SchubertPolyBasis, SepDescPolyBasis
from .quantum_schubert_ring import QDSx, QPDSx, QPSx, QSx, QuantumDoubleSchubertElement, QuantumDoubleSchubertRing, QuantumSingleSchubertRing, make_parabolic_quantum_basis
from .schubert_ring import DoubleSchubertElement, DoubleSchubertRing, DSx, ElemDoubleSchubertRing, SingleSchubertRing, Sx
from .separated_descents import SeparatedDescentsRing, SeparatedDescentsRingElement
from .tensor_ring import TensorRing, TensorRingElement

ASx = FreeAlgebra(basis=SchubertBasis)
ADSx = SeparatedDescentsRing(DSx([]).ring)
J = FreeAlgebra(basis=JBasis)
L = FreeAlgebra(basis=NElementaryBasis)
Z = FreeAlgebra(basis=ZBasis)

__all__ = [
    "FA",
    "PA",
    "ADSx",
    "ASx",
    "DSx",
    "DoubleSchubertElement",
    "DoubleSchubertRing",
    "ElemDoubleSchubertRing",
    "ElemSymPolyBasis",
    "ElementaryBasis",
    "FreeAlgebra",
    "FreeAlgebraBasis",
    "FreeAlgebraElement",
    "J",
    "JBasis",
    "L",
    "MonomialBasis",
    "NElementaryBasis",
    "NilHeckeRing",
    "PolynomialAlgebra",
    "PolynomialAlgebraElement",
    "PolynomialBasis",
    "QDSx",
    "QPDSx",
    "QPSx",
    "QSx",
    "QuantumDoubleSchubertElement",
    "QuantumDoubleSchubertRing",
    "QuantumSingleSchubertRing",
    "SchubertBasis",
    "SchubertPolyBasis",
    "SchubertSchurBasis",
    "SepDescPolyBasis",
    "SeparatedDescentsBasis",
    "SeparatedDescentsRing",
    "SeparatedDescentsRingElement",
    "SingleSchubertRing",
    "Sx",
    "TensorRing",
    "TensorRingElement",
    "WordBasis",
    "Z",
    "ZBasis",
    "make_parabolic_quantum_basis",
]

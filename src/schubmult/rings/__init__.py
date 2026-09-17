"""Ring structures built on combinatorial bases.

Subpackages:

- `schubert`: the Schubert-family rings (``Sx``, ``DSx``, ``Gx``, ``QSx``, ...) -- the
  main user-facing algebra interface.
- `polynomial_algebra`: the polynomial ring ``Z[x_1, x_2, ...]`` with pluggable bases
  (monomial, Schubert, key, slide, Grothendieck, ...).
- `free_algebra`: noncommutative free algebras with combinatorial bases.
- `combinatorial`: rings whose basis elements are combinatorial objects (RC graphs,
  BPDs, plactic classes, ...) rather than permutations.

Modules here: `base_ring` (the shared dict-based ring/element machinery), `tensor_ring`,
`product_ring`, `direct_product_ring`, `printing` (display symbols), `nsym`,
`quasisymmetric_functions`, `thompson_algebra`.

This ``__init__`` re-exports nothing (the commented-out block below is legacy);
import from the subpackages directly.
"""

# from ._separated_descents import SeparatedDescentsRing
# from .free_algebra import FA, FreeAlgebra, FreeAlgebraElement
# from .free_algebra_basis import ElementaryBasis, FreeAlgebraBasis, JBasis, NElementaryBasis, SchubertBasis, SchubertSchurBasis, SeparatedDescentsBasis, WordBasis, ZBasis
# from .nil_hecke import NilHeckeRing
# from .polynomial_algebra import PA, PolynomialAlgebra, PolynomialAlgebraElement
# from .polynomial_basis import ElemSymPolyBasis, MonomialBasis, PolynomialBasis, SchubertPolyBasis, SepDescPolyBasis
# from .schubert.quantum_schubert_ring import QDSx, QPDSx, QPSx, QSx, QuantumDoubleSchubertElement, QuantumDoubleSchubertRing, QuantumSingleSchubertRing, make_parabolic_quantum_basis
# from .schubert.schubert_ring import DoubleSchubertElement, DoubleSchubertRing, DSx, ElemDoubleSchubertRing, SingleSchubertRing, Sx
# from .tensor_ring import TensorRing, TensorRingElement
# from .combinatorial.bpd_ring import BPDRing, BPDRingElement
# from .combinatorial.rc_graph_ring import RCGraphRing, RCGraphRingElement
# from .combinatorial.schubert_monomial_ring import SchubertMonomialPrintingTerm, SchubertMonomialRing, SchubertMonomialRingElement

# ASx = FreeAlgebra(basis=SchubertBasis)
# # Will do ADSx again
# ADSx = SeparatedDescentsRing(DSx([]).ring)
# J = FreeAlgebra(basis=JBasis)
# L = FreeAlgebra(basis=NElementaryBasis)
# Z = FreeAlgebra(basis=ZBasis)

# __all__ = [
#     "FA",
#     "PA",
#     "ADSx",
#     "ASx",
#     "DSx",
#     "DoubleSchubertElement",
#     "DoubleSchubertRing",
#     "ElemDoubleSchubertRing",
#     "ElemSymPolyBasis",
#     "ElementaryBasis",
#     "FreeAlgebra",
#     "FreeAlgebraBasis",
#     "FreeAlgebraElement",
#     "J",
#     "JBasis",
#     "L",
#     "MonomialBasis",
#     "NElementaryBasis",
#     "NilHeckeRing",
#     "PolynomialAlgebra",
#     "PolynomialAlgebraElement",
#     "PolynomialBasis",
#     "QDSx",
#     "QPDSx",
#     "QPSx",
#     "QSx",
#     "QuantumDoubleSchubertElement",
#     "QuantumDoubleSchubertRing",
#     "QuantumSingleSchubertRing",
#     "SchubertBasis",
#     "SchubertPolyBasis",
#     "SchubertSchurBasis",
#     "SepDescPolyBasis",
#     "SeparatedDescentsBasis",
#     "SingleSchubertRing",
#     "Sx",
#     "TensorRing",
#     "TensorRingElement",
#     "WordBasis",
#     "Z",
#     "ZBasis",
#     "make_parabolic_quantum_basis",
# ]

# __all__ = [
#     "BPDRing",
#     "BPDRingElement",
#     "RCGraphRing",
#     "RCGraphRingElement",
#     "SchubertMonomialRing",
#     "SchubertMonomialRingElement",
#     "SchubertMonomialPrintingTerm",
# ]

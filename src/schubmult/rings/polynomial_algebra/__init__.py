"""Polynomial algebra ring and pre-built basis instances.

Exports the core :class:`PolynomialAlgebra` and :class:`PolynomialAlgebraElement`
classes as well as ready-to-use ring instances:

- ``Schub`` — Schubert polynomial basis
- ``Forest`` — Forest polynomial basis
- ``Key`` — Key polynomial (Demazure character) basis
- ``FSlide`` — Fundamental slide polynomial basis
- ``Monomial`` — Standard monomial basis
"""

from ...abc import x, y
from ..schubert.schubert_ring import Sx
from ._core import PA, PolynomialAlgebra, PolynomialAlgebraElement
from .composition_schubert_poly_basis import CompositionSchubertPolyBasis
from .double_forest_poly_basis import DoubleForestPolyBasis
from .forest_poly_basis import ForestPolyBasis
from .key_poly_basis import KeyPolyBasis
from .polynomial_basis import ElemSymPolyBasis, FundamentalSlidePolyBasis, MonomialBasis, MonomialSlidePolyBasis, PolynomialBasis, SchubertPolyBasis, SepDescPolyBasis

Forest = PolynomialAlgebra(ForestPolyBasis(x))
Schub = PolynomialAlgebra(SchubertPolyBasis(x))
Key = PolynomialAlgebra(KeyPolyBasis(x))
FSlide = PolynomialAlgebra(FundamentalSlidePolyBasis(x))
Monomial = PolynomialAlgebra(MonomialBasis(x))
DoubleForest = PolynomialAlgebra(DoubleForestPolyBasis(x, y))

__all__ = [
    "PA",
    "CompositionSchubertPolyBasis",
    "DoubleForest",
    "DoubleForestPolyBasis",
    "ElemSymPolyBasis",
    "FSlide",
    "Forest",
    "ForestPolyBasis",
    "FundamentalSlidePolyBasis",
    "Key",
    "KeyPolyBasis",
    "Monomial",
    "MonomialBasis",
    "MonomialSlidePolyBasis",
    "PolynomialAlgebra",
    "PolynomialAlgebraElement",
    "PolynomialBasis",
    "Schub",
    "SchubertPolyBasis",
    "SepDescPolyBasis",
]

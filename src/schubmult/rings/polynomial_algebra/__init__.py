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
from .glide_poly_basis import GlidePolyBasis
from .grothendieck_poly_basis import GrothendieckPolyBasis
from .grove_poly_basis import GrovePolyBasis
from .key_poly_basis import KeyPolyBasis
from .lascoux_poly_basis import LascouxPolyBasis
from .polynomial_basis import ElemSymPolyBasis, FundamentalSlidePolyBasis, MonomialBasis, MonomialSlidePolyBasis, PolynomialBasis, SchubertPolyBasis, SepDescPolyBasis

ForestPoly = PolynomialAlgebra(ForestPolyBasis(x))
Grove = PolynomialAlgebra(GrovePolyBasis(x))
GrovePoly = PolynomialAlgebra(GrovePolyBasis(x))
Schub = PolynomialAlgebra(SchubertPolyBasis(x))
Key = PolynomialAlgebra(KeyPolyBasis(x))
KeyPoly = PolynomialAlgebra(KeyPolyBasis(x))
FSlide = PolynomialAlgebra(FundamentalSlidePolyBasis(x))
SlidePoly = PolynomialAlgebra(FundamentalSlidePolyBasis(x))
GlidePoly = PolynomialAlgebra(GlidePolyBasis(x))
Monomial = PolynomialAlgebra(MonomialBasis(x))
DoubleForest = PolynomialAlgebra(DoubleForestPolyBasis(x, y))
GrothendieckPoly = PolynomialAlgebra(GrothendieckPolyBasis(x))
LascouxPoly = PolynomialAlgebra(LascouxPolyBasis(x))

__all__ = [
    "PA",
    "CompositionSchubertPolyBasis",
    "DoubleForest",
    "DoubleForestPolyBasis",
    "ElemSymPolyBasis",
    "FSlide",
    "ForestPoly",
    "ForestPolyBasis",
    "FundamentalSlidePolyBasis",
    "GlidePoly",
    "GlidePolyBasis",
    "GrothendieckPoly",
    "GrothendieckPolyBasis",
    "GrovePoly",
    "GrovePolyBasis",
    "Key",
    "KeyPoly",
    "KeyPolyBasis",
    "LascouxPoly",
    "LascouxPolyBasis",
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

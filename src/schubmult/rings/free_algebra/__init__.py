"""The free algebra: the graded dual of `schubmult.rings.polynomial_algebra`.

The free (noncommutative) algebra on generators indexed by nonnegative integers is
the graded dual of the polynomial ring ``Z[x_1, x_2, ...]``, under the pairing in
which the word ``(a_1, ..., a_n)`` is dual to the monomial
``x_1^{a_1} x_2^{a_2} ... x_n^{a_n}`` (i.e. a word is the exponent vector of its dual
monomial). Concretely, ``FreeAlgebraElement.pairing(poly_elem)`` /
``PolynomialAlgebraElement.apply_dual_element(free_elem)`` sum the products of
matching coefficients once both sides are in the word / monomial basis.

Under this duality:

- concatenation of words (the `WordBasis` product) is dual to the variable-splitting
  coproduct on polynomials (`PolynomialAlgebraElement.branch`/``coproduct``), and the
  free algebra's coproduct is dual to polynomial multiplication;
- every free-algebra basis is the dual of a polynomial-algebra basis, exposed via
  ``Basis.dual_basis()``: `SchubertBasis` <-> ``SchubertPolyBasis``, `KeyBasis` <->
  ``KeyPolyBasis``, `ForestBasis` <-> ``ForestPolyBasis``, `FundamentalSlideBasis` <->
  ``FundamentalSlidePolyBasis``, `GrothendieckBasis` <-> ``GrothendieckPolyBasis``, and
  so on, with `WordBasis` <-> ``MonomialBasis`` as the pair everything is computed through.

Basis keys carry a *length* (number of variables) alongside the combinatorial index,
e.g. a `SchubertBasis` key is ``(perm, numvars)``: the dual of ``S_perm`` regarded as
a polynomial in exactly ``numvars`` variables. This grading by number of variables is
what makes the duality with words of a fixed length work.

The core classes are `FreeAlgebra` (the ring, parametrized by a `FreeAlgebraBasis`
class) and `FreeAlgebraElement`; ``change_basis`` converts between bases.

Pre-built instances:
    - ``FA``: `WordBasis` (the default)
    - ``ASx``: `SchubertBasis`
    - ``AGx``: `GrothendieckBasis`
    - ``ADSx``: the double Schubert separated-descents ring used for expansion
    - ``ForestDual``, ``GroveDual``, ``GlideDual``: the corresponding bases

Available bases:
    WordBasis, SchubertBasis, CompositionSchubertBasis, ElementaryBasis,
    ForestBasis, FundamentalSlideBasis, GlideBasis, GrothendieckBasis, GroveBasis,
    JBasis, JTBasis, KeyBasis, LascouxBasis, MonomialSlideBasis, NElementaryBasis,
    SchubertSchurBasis, SchurElementaryBasis, SeparatedDescentsBasis, ZBasis.
"""

from ._core import FA, ADSx, AGx, ASx, FreeAlgebra, FreeAlgebraElement, make_AGx
from .composition_schubert_basis import CompositionSchubertBasis
from .elementary_basis import ElementaryBasis
from .forest_basis import ForestBasis, ForestDual
from .free_algebra_basis import FreeAlgebraBasis
from .fundamental_slide_basis import FundamentalSlideBasis
from .glide_basis import GlideBasis
from .grothendieck_basis import GrothendieckBasis
from .grove_basis import GroveBasis
from .j_basis import JBasis
from .jt_basis import JTBasis
from .key_basis import KeyBasis
from .lascoux_basis import LascouxBasis
from .monomial_slide_basis import MonomialSlideBasis
from .nelementary_basis import NElementaryBasis
from .schubert_basis import SchubertBasis
from .schubert_schur_basis import SchubertSchurBasis
from .schur_elementary_basis import SchurElementaryBasis
from .separated_descents_basis import SeparatedDescentsBasis, _SeparatedDescentsBasis
from .word_basis import WordBasis
from .z_basis import ZBasis

GroveDual = FreeAlgebra(GroveBasis)
GlideDual = FreeAlgebra(GlideBasis)

__all__ = [
    "FA",
    "ADSx",
    "AGx",
    "ASx",
    "CompositionSchubertBasis",
    "ElementaryBasis",
    "ForestBasis",
    "ForestDual",
    "FreeAlgebra",
    "FreeAlgebraBasis",
    "FreeAlgebraElement",
    "FundamentalSlideBasis",
    "GlideBasis",
    "GlideDual",
    "GrothendieckBasis",
    "GroveBasis",
    "GroveDual",
    "JBasis",
    "JTBasis",
    "KeyBasis",
    "LascouxBasis",
    "MonomialSlideBasis",
    "NElementaryBasis",
    "SchubertBasis",
    "SchubertSchurBasis",
    "SchurElementaryBasis",
    "SeparatedDescentsBasis",
    "WordBasis",
    "ZBasis",
    "_SeparatedDescentsBasis",
    "make_AGx",
]

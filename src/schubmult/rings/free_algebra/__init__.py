"""Free algebra module providing multiple bases for Schubert calculus.

The core classes are :class:`FreeAlgebra` (the ring) and
:class:`FreeAlgebraElement` (its elements).  Elements can be expressed in
any of the available bases and converted between them via ``change_basis``.

Pre-built instances:
    - ``FA``: FreeAlgebra with WordBasis (default)
    - ``ASx``: FreeAlgebra with SchubertBasis
    - ``ADSx``: FreeAlgebra with double Schubert basis

Available bases:
    WordBasis, SchubertBasis, CompositionSchubertBasis, ElementaryBasis,
    ForestBasis, FundamentalSlideBasis, JBasis, JTBasis, KeyBasis,
    MonomialSlideBasis, NElementaryBasis, SchubertSchurBasis,
    SeparatedDescentsBasis, ZBasis.
"""

from ._core import FA, ADSx, ASx, FreeAlgebra, FreeAlgebraElement
from .composition_schubert_basis import CompositionSchubertBasis
from .elementary_basis import ElementaryBasis
from .forest_basis import ForestBasis, ForestDual
from .free_algebra_basis import FreeAlgebraBasis
from .fundamental_slide_basis import FundamentalSlideBasis
from .j_basis import JBasis
from .jt_basis import JTBasis
from .key_basis import KeyBasis
from .monomial_slide_basis import MonomialSlideBasis
from .nelementary_basis import NElementaryBasis
from .schubert_basis import SchubertBasis
from .schubert_schur_basis import SchubertSchurBasis
from .separated_descents_basis import SeparatedDescentsBasis, _SeparatedDescentsBasis
from .word_basis import WordBasis
from .z_basis import ZBasis

__all__ = [
    "FA",
    "ADSx",
    "ASx",
    "CompositionSchubertBasis",
    "ElementaryBasis",
    "ForestBasis",
    "ForestDual",
    "FreeAlgebra",
    "FreeAlgebraBasis",
    "FreeAlgebraElement",
    "FundamentalSlideBasis",
    "JBasis",
    "JTBasis",
    "KeyBasis",
    "MonomialSlideBasis",
    "NElementaryBasis",
    "SchubertBasis",
    "SchubertSchurBasis",
    "SeparatedDescentsBasis",
    "WordBasis",
    "ZBasis",
    "_SeparatedDescentsBasis",
]

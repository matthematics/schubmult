from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import S
from schubmult.utils.tuple_utils import pad_tuple

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis
from .schubert_basis import SchubertBasis


class CompositionSchubertBasis(FreeAlgebraBasis):
    """Schubert basis indexed by padded trimcode compositions.

    A key is a composition ``c`` and corresponds to Schubert key
    ``(uncode(c), len(c))``.
    """

    @classmethod
    def is_key(cls, x):
        return isinstance(x, (tuple, list))

    @classmethod
    def as_schubert_key(cls, key):
        return (uncode(key), len(key))

    @classmethod
    def as_key(cls, key):
        if isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], int) and isinstance(key[0], (Permutation, list, tuple)):
            perm = Permutation(key[0])
            return pad_tuple(perm.trimcode, key[1])
        return tuple(key)

    @classmethod
    def from_rc_graph(cls, rc_graph):
        out = {}
        for key, coeff in SchubertBasis.from_rc_graph(rc_graph).items():
            out[cls.as_key(key)] = coeff
        return out

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        return {cls.as_key(k): v for k, v in SchubertBasis.product(cls.as_schubert_key(key1), cls.as_schubert_key(key2), coeff).items()}

    zero_monom = ()

    @classmethod
    def coproduct(cls, key):
        return {cls.as_key(k): v for k, v in SchubertBasis.coproduct(cls.as_schubert_key(key)).items()}

    @classmethod
    def bcoproduct(cls, key):
        return {cls.as_key(k): v for k, v in SchubertBasis.bcoproduct(cls.as_schubert_key(key)).items()}

    @classmethod
    def internal_product(cls, key1, key2, coeff=S.One):
        return {cls.as_key(k): v for k, v in SchubertBasis.internal_product(cls.as_schubert_key(key1), cls.as_schubert_key(key2), coeff).items()}

    @classmethod
    def skew_element(cls, w, u, n):
        return {cls.as_key(k): v for k, v in SchubertBasis.skew_element(w, u, n).items()}


    @classmethod
    def dual_basis(cls):
        return SchubertBasis.dual_basis()

    @classmethod
    def transition(cls, other_basis):
        if other_basis == CompositionSchubertBasis:
            return lambda x: {x: S.One}
        if other_basis == SchubertBasis:
            return lambda x: {cls.as_schubert_key(x): S.One}
        return lambda x: SchubertBasis.transition(other_basis)(cls.as_schubert_key(x))

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(cls.as_key(k), "CompSchub")

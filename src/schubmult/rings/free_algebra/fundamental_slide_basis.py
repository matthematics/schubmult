from functools import cache

from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class FundamentalSlideBasis(FreeAlgebraBasis):
    @staticmethod
    def _weak_compositions(length, total):
        if length == 0:
            if total == 0:
                yield ()
            return
        if length == 1:
            yield (total,)
            return
        for i in range(total + 1):
            for tail in FundamentalSlideBasis._weak_compositions(length - 1, total - i):
                yield (i, *tail)

    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        return tuple(x)

    # @classmethod
    # def from_rc_graph(cls, rc_graph):
    #     return {rc_graph.length_vector(): 1}

    # @classmethod
    # def product(cls, key1, key2, coeff=S.One):
    #     return {(*key1, *key2): coeff}

    zero_monom = ()


    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "FS")


    @classmethod
    def transition_schubert(cls, key):
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(WordBasis.transition_schubert, cls.transition_word(key))

    @classmethod
    @cache
    def transition_word(cls, key):
        from schubmult.abc import x
        from schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis import FundamentalSlidePolyBasis
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        key = tuple(key)
        length = len(key)
        total = sum(key)
        mon = MonomialBasis(x)
        fslide = FundamentalSlidePolyBasis(x)

        ret = {}
        for word_key in cls._weak_compositions(length, total):
            dct = mon.transition_slide({word_key: S.One}, fslide)
            coeff = dct.get(key, S.Zero)
            if coeff != S.Zero:
                ret[word_key] = coeff
        return ret


    @classmethod
    def dual_basis(cls):
        from ..polynomial_algebra.fundamental_slide_poly_basis import FundamentalSlidePolyBasis
        return FundamentalSlidePolyBasis

    @classmethod
    def transition(cls, other_basis):
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))

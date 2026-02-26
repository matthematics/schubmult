from functools import cache
from itertools import combinations

from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class MonomialSlideBasis(FreeAlgebraBasis):
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
            for tail in MonomialSlideBasis._weak_compositions(length - 1, total - i):
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
        return GenericPrintingTerm(k, "MS")


    @classmethod
    @cache
    def transition_fundamental_slide(cls, key):
        key = tuple(key)
        n = len(key)
        flat_key = tuple(c for c in key if c > 0)

        if len(flat_key) == 0:
            return {tuple([0] * n): S.One}

        prefix_key = []
        running = 0
        for value in key:
            running += value
            prefix_key.append(running)

        def coarsenings(comp):
            if len(comp) == 1:
                return [comp]

            ret = []

            def recurse(index, current_sum, current_parts):
                if index == len(comp):
                    ret.append((*current_parts, current_sum))
                    return
                recurse(index + 1, current_sum + comp[index], current_parts)
                recurse(index + 1, comp[index], (*current_parts, current_sum))

            recurse(1, comp[0], ())
            return ret

        ret = {}

        for coarse in coarsenings(flat_key):
            m = len(coarse)
            for positions in combinations(range(n), m):
                weak_comp = [0] * n
                for i, pos in enumerate(positions):
                    weak_comp[pos] = coarse[i]

                running = 0
                valid = True
                for i, value in enumerate(weak_comp):
                    running += value
                    if running > prefix_key[i]:
                        valid = False
                        break

                if valid:
                    ret[tuple(weak_comp)] = S.One

        return ret


    @classmethod
    def transition(cls, other_basis):
        from .fundamental_slide_basis import FundamentalSlideBasis
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        if other_basis == WordBasis:
            return lambda x: cls.transition_word(x)

        if other_basis == FundamentalSlideBasis:
            return lambda x: cls.transition_fundamental_slide(x)
        if other_basis == SchubertBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition_schubert, cls.transition_word(x))
        return lambda x: FreeAlgebraBasis.compose_transition(FundamentalSlideBasis.transition(other_basis), cls.transition_fundamental_slide(x))

    @classmethod
    @cache
    def transition_word(cls, key):
        from schubmult.abc import x
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis
        from schubmult.rings.polynomial_algebra.monomial_slide_poly_basis import MonomialSlidePolyBasis

        key = tuple(key)
        length = len(key)
        total = sum(key)
        mon = MonomialBasis(x)
        mslide = MonomialSlidePolyBasis(x)

        ret = {}
        for word_key in cls._weak_compositions(length, total):
            dct = mon.transition_slide({word_key: S.One}, mslide)
            coeff = dct.get(key, S.Zero)
            if coeff != S.Zero:
                ret[word_key] = coeff
        return ret

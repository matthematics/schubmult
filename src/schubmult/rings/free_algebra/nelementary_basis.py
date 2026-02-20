from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class NElementaryBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        return tuple(x)

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        return {(*key1, *key2): coeff}

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "L")

    @classmethod
    def transition_word(cls, tup):
        from sage.all import Composition

        pain = Composition(tup)
        piss = list(pain.finer())
        return {tuple([int(p) for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

    @classmethod
    def transition(cls, other_basis):
        from .word_basis import WordBasis

        return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(x))

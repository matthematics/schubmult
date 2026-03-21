from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class NElementaryBasis(FreeAlgebraBasis):
    """Non-commutative elementary basis (L basis) of the free algebra.

    Keys are tuples of positive integers. The transition to the word
    basis uses composition refinements from SageMath.
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a tuple or list."""
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* to a tuple key."""
        return tuple(x)

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        """Concatenate two keys."""
        return {(*key1, *key2): coeff}

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        """Return an ``L``-labelled display object for key *k*."""
        return GenericPrintingTerm(k, "L")

    @classmethod
    def transition_word(cls, tup):
        """Transition an NElementary key to the word basis via composition refinements (requires SageMath)."""
        from sage.all import Composition

        pain = Composition(tup)
        piss = list(pain.finer())
        return {tuple([int(p) for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from NElementaryBasis to *other_basis*."""
        from .word_basis import WordBasis

        return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(x))

from functools import cache

from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import Symbol, sstr

from ..schubert.schubert_ring import Sx
from .free_algebra_basis import FreeAlgebraBasis


class SchubertSchurBasis(FreeAlgebraBasis):
    """Schubert-Schur basis of the free algebra.

    Keys are ``(partition_tuple, Permutation)`` pairs encoding products of
    a Grassmannian Schubert polynomial with a Schur polynomial.
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a ``(list/tuple, Permutation/list/tuple)`` pair."""
        return len(x) == 2 and isinstance(x[0], list | tuple) and isinstance(x[1], Permutation | list | tuple)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* into a ``(tuple, Permutation)`` key."""
        return (tuple(x[0]), Permutation(x[1]))

    zero_monom = ((), Permutation([]))

    @classmethod
    @cache
    def coproduct(cls, key):
        """Compute the coproduct of a Schubert-Schur key via the Schubert basis."""
        from ...utils._mul_utils import _tensor_product_of_dicts_first
        from .schubert_basis import SchubertBasis

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(SchubertBasis.transition(cls)(x[0]), SchubertBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda x: SchubertBasis.coproduct(x), cls.transition_schubert(*key)),
        )

    @classmethod
    def transition_schubert(cls, lambd, perm):
        """Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis."""
        from schubmult.symbolic.poly.variables import MaskedGeneratingSet

        from ..schubert.schubert_ring import SingleSchubertRing

        if lambd[-1] == 0:
            return {(perm, len(lambd)): 1}
        numvars = len(lambd)
        extra = lambd[-1] + len(lambd) - 1
        dom = uncode([numvars] * extra)
        grass_perm = uncode(lambd) * dom
        w0 = uncode(list(range(numvars - 1, 0, -1)))
        lower_perm = perm * w0
        dom_perm = uncode(([numvars] * extra) + list(range(numvars - 1, 0, -1)))
        shifted_ring = SingleSchubertRing(MaskedGeneratingSet(Sx([]).ring.genset, list(range(1, extra + 1))))
        start_schub = Sx(grass_perm)
        start_schub *= shifted_ring(lower_perm).in_SEM_basis()
        return {(k * ~dom_perm, numvars): v for k, v in start_schub.items()}

    @classmethod
    def transition_word(cls, lambd, perm):
        """Transition a Schubert-Schur key to the word basis via the Schubert basis."""
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis), cls.transition_schubert(lambd, perm))

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from SchubertSchurBasis to *other_basis*."""
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: x
        return None

    @classmethod
    def printing_term(cls, k):
        """Return an ``SS``-prefixed symbol for key *k*."""
        return Symbol(f"SS{sstr(k)}")

from functools import cache

from schubmult.combinatorics.permutation import uncode

from ..printing import GenericPrintingTerm
from ..schubert.schubert_ring import Sx
from .free_algebra_basis import FreeAlgebraBasis


class SchurElementaryBasis(FreeAlgebraBasis):
    """Schur-Elementary basis of the free algebra.

    Keys are ``(tuple, tuple)`` pairs encoding products of
    a standard elementary monomial (first tuple) and a Schur polynomial (second tuple, partition
    of length percisely len(first_tuple) + 1 in increasing order, with zeros at the beginning if
    needed).  For example, the key ``((1, 2), (0, 1, 3))`` corresponds to the product of the elementary monomial
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a ``(list/tuple, list/tuple)`` pair."""
        return len(x) == 2 and isinstance(x[0], list | tuple) and isinstance(x[1], list | tuple)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* into a ``(tuple, tuple)`` key."""
        return (tuple(x[0]), tuple(x[1]))

    zero_monom = ((), ())

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
    def transition_schubert(cls, elem_tup, lambd):
        """Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis."""
        from schubmult.symbolic import expand_seq

        from .elementary_basis import ElementaryBasis

        numvars = len(lambd)
        if lambd[-1] == 0:
            return ElementaryBasis.transition_schubert((*elem_tup, 0), numvars)
        extra = lambd[-1] + len(lambd) - 1
        dom = uncode([numvars] * extra)
        grass_perm = uncode(lambd) * dom
        dom_perm = uncode(([numvars] * extra) + list(range(numvars - 1, 0, -1)))
        #shifted_ring = SingleSchubertRing(MaskedGeneratingSet(Sx([]).ring.genset, list(range(1, extra + 1))))
        start_schub = Sx(grass_perm)
        start_schub *= Sx.from_expr(expand_seq([numvars - 1 - i - a for i, a in enumerate(reversed(elem_tup))], Sx.genset[extra+1:]))
        return {(k * ~dom_perm, numvars): v for k, v in start_schub.items()}

    @classmethod
    def transition_word(cls, elem_tup, lambd):
        """Transition a Schubert-Schur key to the word basis via the Schubert basis."""
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis), cls.transition_schubert(elem_tup, lambd))

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from SchurElementaryBasis to *other_basis*."""
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if other_basis == SchurElementaryBasis:
            return lambda x: x
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(*x))

    @classmethod
    def printing_term(cls, k):
        """Return an ``SE``-prefixed symbol for key *k*."""
        #return Symbol(f"SE{sstr(k)}")
        return GenericPrintingTerm(k, "SE")

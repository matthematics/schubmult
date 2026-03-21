from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from ..schubert.schubert_ring import Sx
from .free_algebra_basis import FreeAlgebraBasis


class ElementaryBasis(FreeAlgebraBasis):
    """Elementary symmetric function basis of the free algebra.

    Keys are ``(tuple, int)`` pairs where the tuple encodes an elementary
    symmetric function composition and the integer is the number of variables.
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a ``(tuple/list, int)`` pair."""
        return isinstance(x, tuple | list) and len(x) == 2 and isinstance(x[0], tuple | list) and isinstance(x[1], int)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* into a ``(tuple, int)`` key."""
        return (x[0], x[1])

    zero_monom = ((), 0)

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from ElementaryBasis to *other_basis*."""
        from .schubert_basis import SchubertBasis
        from .schubert_schur_basis import SchubertSchurBasis
        from .word_basis import WordBasis

        if other_basis == cls:
            return lambda x: {x: S.One}
        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_schubert_schur(*y), cls.transition_schubert(*x))
        if other_basis == WordBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_word(*y), cls.transition_schubert(*x))
        if other_basis.__name__ == "_SeparatedDescentsBasis":
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_separated_descents(other_basis.k, *y), cls.transition_schubert(*x))
        return None

    @classmethod
    def transition_schubert(cls, tup, numvars):
        """Transition an elementary key to the Schubert basis."""
        from schubmult.abc import x
        from schubmult.symbolic import prod
        from schubmult.symbolic.poly.schub_poly import monom_sym

        mu = list(range(numvars, 0, -1))
        if len(mu) < len(tup):
            mu = [*([numvars] * (len(tup) - len(mu))), *mu]
        tup = tuple(reversed(tup))
        flat_part = tup[: -numvars + 1]
        boink_part = tup[-numvars + 1 :]
        painted_bagel = Sx([]).ring.zero
        pickles = [mu[i] - flat_part[len(flat_part) - 1 - i] for i in range(len(flat_part))]
        painted_bagel = Sx.from_expr(monom_sym(pickles, len(flat_part), Sx([]).ring.genset))

        painted_bagel *= prod([x[i + 1] ** (mu[i] - boink_part[i - len(flat_part)]) for i in range(len(flat_part), len(mu))])
        w0 = ~uncode(mu)
        monom = {}
        for k, v in painted_bagel.items():
            if (k * w0).inv != w0.inv - k.inv:
                raise Exception
            monom[(k * w0, numvars)] = v
        return dict(monom)

    @classmethod
    def printing_term(cls, k):
        """Return an ``Elem``-labelled display object for key *k*."""
        return GenericPrintingTerm(k, "Elem")

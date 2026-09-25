"""`ElementaryBasis`: free-algebra basis indexed by ``(composition, numvars)``, dual to products of
elementary symmetric polynomials ``e_{c_1}(x_1..x_k) e_{c_2}(x_1..x_{k-1}) ...`` in nested
variable sets. `SchubertBasis` expands into it via the monomials of ``S_{perm * w0}``.
"""

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
        if cls.is_key(x):
            return (tuple(x[0]), x[1])
        return None

    zero_monom = ((), 0)

    @classmethod
    def canonical_key(cls, tup, numvars):
        """Canonicalize ``(tup, numvars)``: the flag part ``tup[:numvars-1]`` is kept as is; the
        symmetric tail ``tup[numvars-1:]`` (indices of ``e_k(x_1..x_numvars)`` factors) is sorted
        with zeros dropped, or ``(0,)`` if empty.
        """
        if numvars == 0:
            return ((), 0)
        head = tuple(tup[: numvars - 1])
        tail = tuple(sorted(t for t in tup[numvars - 1 :] if t))
        return ((*head, *tail) if tail else (*head, 0), numvars)

    @staticmethod
    def staircase(numvars, degree):
        """Dominant code ``[n]*(L-1) + [n, n-1, ..., 1]`` with ``L = max(degree, 1)`` symmetric slots.

        Every degree-``degree`` elementary product in ``numvars`` variables has at most ``degree``
        tail factors, so using ``L`` slots uniformly makes the Cauchy-kernel duality cover the
        whole degree at once; a per-key tail length only gives duality within a smaller span.
        """
        L = max(degree, 1)
        return L, [numvars] * (L - 1) + list(range(numvars, 0, -1))

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
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(*x))

    @classmethod
    def transition_schubert(cls, tup, numvars):
        """Transition an elementary key to the Schubert basis."""
        from schubmult.abc import x
        from schubmult.symbolic import prod
        from schubmult.symbolic.common_polys import monom_sym

        if numvars == 0:
            return {(uncode([]), 0): S.One} if not any(tup) else {}
        head = list(tup[: numvars - 1])
        tail = sorted(t for t in tup[numvars - 1 :] if t)
        L, mu = cls.staircase(numvars, sum(tup))
        tail = [0] * (L - len(tail)) + tail
        # symmetric tail lives in the first L variables, flag part in the remaining numvars-1
        painted_bagel = monom_sym([numvars - t for t in tail], L, Sx([]).ring.genset)
        painted_bagel *= prod([x[L + j + 1] ** ((numvars - 1 - j) - head[numvars - 2 - j]) for j in range(numvars - 1)])
        painted_bagel = Sx.from_expr(painted_bagel)
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

    @classmethod
    def dual_basis(cls):
        from ..polynomial_algebra.elem_sym_poly_basis import ElemSymPolyBasis

        return ElemSymPolyBasis

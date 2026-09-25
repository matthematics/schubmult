"""`ElementaryBasis`: free-algebra basis indexed by ``(composition, numvars)``, dual to products of
elementary symmetric polynomials ``e_{a_1}(x_1) e_{a_2}(x_1, x_2) ... e_{a_{n-1}}(x_1..x_{n-1})`` times a
symmetric tail of ``e_k(x_1..x_n)`` factors (`ElemSymPolyBasis`).

Transitions to and from `SchubertBasis` go through the finite ``(numvars, degree)`` block: expanding
every elementary product of that degree in the Schubert basis (Pieri rule) gives the Schubert -> Elem
matrix, and its inverse is Elem -> Schubert.
"""

import itertools
from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from ..schubert.schubert_ring import Sx
from .free_algebra_basis import FreeAlgebraBasis


def _partitions(total, max_part):
    """Partitions of ``total`` with parts in ``1..max_part``, as ascending tuples."""
    if total == 0:
        yield ()
        return
    for part in range(1, min(total, max_part) + 1):
        for rest in _partitions(total - part, part):
            yield (*rest, part)


def _compositions(total, length):
    """Weak compositions of ``total`` into ``length`` parts."""
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, length - 1):
            yield (first, *rest)


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
    def degree_keys(numvars, degree):
        """All canonical keys of the given degree in ``numvars`` variables: flag part ``a_i <= i``,
        tail a partition with parts in ``1..numvars``. There are as many as monomials of that degree.
        """
        if numvars == 0:
            return [((), 0)] if degree == 0 else []
        keys = []
        for head in itertools.product(*[range(i + 1) for i in range(1, numvars)]):
            rest = degree - sum(head)
            if rest < 0:
                continue
            if rest == 0:
                keys.append(((*head, 0), numvars))
            else:
                keys.extend(((*head, *tail), numvars) for tail in _partitions(rest, numvars))
        return keys

    @staticmethod
    def _elem_product_schubert(key):
        """Schubert expansion of the elementary product ``E_key`` (Pieri rule)."""
        from schubmult.abc import e

        tup, numvars = key
        res = Sx.one
        for i, a in enumerate(tup):
            if a:
                res = Sx.elem_mul(res, e(a, min(i + 1, numvars), Sx.genset[1:]))
        return res

    @classmethod
    @cache
    def schubert_block(cls, numvars, degree):
        """``(keys, perms, to_schubert, to_elementary)`` for one ``(numvars, degree)`` block.

        ``perms`` are the permutations whose Schubert polynomial lies in ``Z[x_1..x_numvars]`` with
        that degree (Lehmer codes of length ``numvars``). ``to_schubert[key][perm]`` is the
        coefficient of ``S_perm`` in ``E_key``; ``to_elementary[perm][key]`` is the inverse matrix,
        i.e. the coefficient of ``Elem(key)`` in ``Schub(perm)``.
        """
        from sympy import Matrix

        keys = cls.degree_keys(numvars, degree)
        perms = [uncode(list(code)) for code in _compositions(degree, numvars)] if numvars else [uncode([])]
        perm_index = {perm: i for i, perm in enumerate(perms)}
        to_schubert = {}
        rows = [[0] * len(keys) for _ in perms]
        for j, key in enumerate(keys):
            expansion = {perm: int(c) for perm, c in cls._elem_product_schubert(key).items() if c != 0}
            to_schubert[key] = expansion
            for perm, c in expansion.items():
                rows[perm_index[perm]][j] = c
        inverse = Matrix(rows).inv()
        to_elementary = {perm: {keys[j]: int(inverse[j, i]) for j in range(len(keys)) if inverse[j, i] != 0} for i, perm in enumerate(perms)}
        return keys, perms, to_schubert, to_elementary

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
        """Transition an elementary key to the Schubert basis (row of the inverse block matrix)."""
        if numvars == 0:
            return {(uncode([]), 0): S.One} if not any(tup) else {}
        key = cls.canonical_key(tup, numvars)
        _keys, perms, _to_schubert, to_elementary = cls.schubert_block(numvars, sum(tup))
        return {(perm, numvars): S(c) for perm in perms if (c := to_elementary[perm].get(key, 0)) != 0}

    @classmethod
    def transition_from_schubert(cls, perm, numvars):
        """Expand ``Schub(perm, numvars)`` in this basis: the Schubert coefficients of each ``E_key``."""
        if numvars == 0:
            return {((), 0): S.One} if perm.inv == 0 else {}
        if len(perm.trimcode) > numvars:
            raise ValueError(f"S_{perm} is not a polynomial in {numvars} variables")
        keys, _perms, to_schubert, _to_elementary = cls.schubert_block(numvars, perm.inv)
        return {key: S(c) for key in keys if (c := to_schubert[key].get(perm, 0)) != 0}

    @classmethod
    def printing_term(cls, k):
        """Return an ``Elem``-labelled display object for key *k*."""
        return GenericPrintingTerm(k, "Elem")

    @classmethod
    def dual_basis(cls):
        from ..polynomial_algebra.elem_sym_poly_basis import ElemSymPolyBasis

        return ElemSymPolyBasis

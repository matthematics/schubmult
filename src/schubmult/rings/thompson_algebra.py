"""`ThompsonAlgebra`: a noncommutative algebra on words in generators ``T_i`` and ``R_i``.

A monomial is a tuple of nonzero integers: ``i > 0`` stands for ``T_i`` and ``i < 0`` for
``R_{-i}``. Products are rewritten to a normal form by the rules in `ThompsonAlgebra._commute_pair`:
``T_i T_j = T_j T_{i+1}`` for ``i > j`` (the Thompson monoid relation), and ``T_i R_j`` moves
``R`` to the left with the index shifts (and one two-term case) given there.
"""

from functools import cache

from sympy import Symbol, UnevaluatedExpr, sympify

from schubmult.symbolic import CoercionFailed, prod
from schubmult.utils._mul_utils import add_perm_dict

from .base_ring import BaseRing, BaseRingElement


class ThompsonAlgebraElement(BaseRingElement):
    """Element of `ThompsonAlgebra`: a dict from normal-form words to coefficients."""

    # def __repr__(self):
    #     return f"ThompsonAlgebraElement({self.data})"

    # def __str__(self):
    #     return str(self.data)

class ThompsonAlgebra(BaseRing):
    """Algebra on words in ``T_i`` (positive index) and ``R_i`` (negative index). See the module docstring."""

    _t = tuple([Symbol(f"T_{i}", commutative=False) for i in range(100)])
    _r = tuple([Symbol(f"R_{i}", commutative=False) for i in range(100)])

    def __init__(self):
        super().__init__()
        self.zero_monom = ()
        self.dtype = type("ThompsonAlgebraElement", (ThompsonAlgebraElement,), {"ring": self})
        #self.zero = self.dtype()
        #self.one = self.from_dict({self.zero_monom: 1})

    def __hash__(self):
        return hash((self.__class__.__name__, tuple(self._t), tuple(self._r)))

    @cache
    def printing_term(self, monomial):
        """Noncommutative product of the ``T_i``/``R_i`` symbols for the word."""
        return UnevaluatedExpr(prod(self._t[i] if i > 0 else self._r[-i] for i in monomial))

    def _commute_pair(self, i, j):
        """Rewrite the adjacent pair ``(i, j)`` (``i`` on the left) into a list of replacement pairs.

        ``T_i T_j -> T_j T_{i+1}`` when ``i > j > 0``; ``T_i R_j`` becomes ``R_{j-1} T_i`` if ``i < j - 1``,
        ``R_{i+1} T_i + R_i T_{i+1}`` if ``i == j - 1``, and ``R_j T_{i+1}`` otherwise. Anything else is
        already in normal form.
        """
        if (i > 0 and j > 0) and i > j:
            return [(j, i + 1)]
        if j < 0 and i > 0:
            real_j = -j
            if i < real_j - 1:
                return [(-(real_j - 1), i)]
            if i == real_j - 1:
                return [(-(i + 1), i), (-i, i + 1)]
            return [(-real_j, i + 1)]
        return [(i, j)]

    def new(self, x):
        """Build an element from a word (normalized via `_mul_monomials`), a number, or an existing element."""
        if isinstance(x, list | tuple):
            return self.from_dict(self._mul_monomials(tuple(x), ()))
            #from_dict({tuple(x): 1})
        if sympify(x).is_Number:
            return self.from_dict({self.zero_monom: x})
        if isinstance(x, ThompsonAlgebraElement) and x.ring == self:
            return x
        raise ValueError(f"Cannot create ThompsonAlgebraElement from {x} of type {type(x)}")


    def _mul_monomials(self, mon1, mon2):
        """Normal form of the concatenation ``mon1 mon2``: peel letters off the end of ``mon1`` and push
        each through ``mon2`` with `_commute_pair`, branching on two-term rewrites.
        """
        if len(mon1) == 0:
            return {tuple(mon2): 1}
        last_elem = mon1[-1]
        ret_mon1 = mon1[:-1]
        if len(mon2) == 0:
            return self._mul_monomials(ret_mon1, [last_elem])
        first_elem = mon2[0]
        commute_result = self._commute_pair(last_elem, first_elem)
        if len(commute_result) == 1:
            ret_mon2 = [*commute_result[0], *mon2[1:]]
            return self._mul_monomials(ret_mon1, ret_mon2)
        ret_mon2_first = [*commute_result[0], *mon2[1:]]
        ret_mon2_second = [*commute_result[1], *mon2[1:]]
        return add_perm_dict(self._mul_monomials(ret_mon1, ret_mon2_first),
                             self._mul_monomials(ret_mon1, ret_mon2_second))

    def mul(self, elem, other):
        """Scalar multiplication, or the bilinear extension of `_mul_monomials`."""
        try:
            other = self.domain_new(other)
            return self.from_dict({k: other * v for k, v in elem.items()})
        except CoercionFailed:
            ret_elem = self.zero
            for mon1, coeff in elem.items():
                for mon2, coeff2 in other.items():
                    new_mon = self._mul_monomials(mon1, mon2)
                    new_coeff = coeff * coeff2
                    ret_elem += new_coeff * self.from_dict(new_mon)
            return ret_elem

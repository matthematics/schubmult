"""Parabolic quantum (single) Schubert polynomial ring: the ``QPSx`` interface.

`ParabolicQuantumSingleSchubertRing` is a `ParabolicQuantumDoubleSchubertRing`
with a zero coefficient alphabet, indexed by a composition ``index_comp``
specifying the parabolic subgroup's block sizes.
"""

from bisect import bisect_left
from functools import cache

import schubmult.mult.quantum as py
import schubmult.mult.quantum_double as yz
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S
from schubmult.symbolic.common_polys import elem_sym_poly, schubpoly_from_elems, xreplace_genvars
from schubmult.symbolic.poly.variables import GeneratingSet, poly_genset
from schubmult.utils.perm_utils import is_parabolic

from .base_schubert_ring import BaseSchubertElement
from .parabolic_quantum_double_schubert_ring import ParabolicQuantumDoubleSchubertElement, ParabolicQuantumDoubleSchubertRing, q_var


class ParabolicQuantumSingleSchubertRing(ParabolicQuantumDoubleSchubertRing):
    """Quantum Schubert polynomials for the partial flag variety with block sizes ``index_comp``.
    Construct via ``QPSx(*index_comp)``; basis permutations must be parabolic for the given blocks.
    """

    def __init__(self, genset, index_comp):
        super().__init__(genset, poly_genset(0), index_comp)

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, self.index_comp, "PQSB"))

    @cache
    def cached_schubpoly(self, k):
        """The explicit parabolic quantum Schubert polynomial for ``k``, expanded against the appropriate
        longest element (extended if ``k`` is larger than the ring's default).
        """
        if len(k) > len(self._longest):
            parabolic_index = []
            start = 0
            index_comp = [*self._n, len(k) + 1 - self._N[-1]]
            for i in range(len(index_comp)):
                end = start + index_comp[i]
                parabolic_index += list(range(start + 1, end))
                start = end
            otherlong = Permutation(list(range(parabolic_index[-1] + 1, 0, -1)))
            longpar = Permutation.longest_element(*parabolic_index)
            longest = otherlong * longpar
        else:
            longest = self._longest
        return schubpoly_from_elems(k, self.genset, self.coeff_genset, elem_func=self.elem_sym, mumu=~longest)

    def elem_sym(self, p, k, varl1, varl2):
        """Parabolic quantum elementary symmetric polynomial ``E_p(x_1..x_k)``: classical below the first
        block boundary, with a ``q``-correction term at each block boundary ``N_j``.
        """
        if p < 0 or p > k:
            return 0
        if p == 0 and k >= 0:
            return 1
        if k <= self._N[1]:
            return elem_sym_poly(p, k, varl1, poly_genset(0))
        ret = 0
        j = bisect_left(self._N, k)
        if j < len(self._N) and k == self._N[j]:
            ret = (-((-1) ** (self._n[j - 1]))) * q_var[j - 1] * self.elem_sym(p - self._N[j] + self._N[j - 2], self._N[j - 2], varl1, varl2)
        ret += self.elem_sym(p, k - 1, varl1, varl2) + varl1[k - 1] * self.elem_sym(p - 1, k - 1, varl1, varl2)
        return ret

    def _coerce_mul(self, other):
        """Accept elements of this ring, or of the double parabolic ring with matching blocks (coefficient alphabet zeroed)."""
        if isinstance(other, BaseSchubertElement):
            if other.ring == self:
                if self.genset == other.ring.genset:
                    return other
            if type(other.ring) is ParabolicQuantumDoubleSchubertRing:
                if self.genset == other.ring.genset and self.index_comp == other.ring.index_comp:
                    newbasis = ParabolicQuantumDoubleSchubertRing(self.genset, poly_genset(0), self.index_comp)
                    return newbasis.from_dict(other)
        return None

    @property
    def coeff_genset(self):
        """Always the zero alphabet."""
        return poly_genset(0)

    @cache
    def cached_product(self, u, v, basis2):
        """Full-flag quantum product (``schubmult_q_fast`` or the generic double kernel), then projected
        to the parabolic ring via ``process_coeff_dict`` (Peterson-Woodward).
        """
        if self == basis2:
            initial_dict = py.schubmult_q_fast({u: S.One}, v)
        else:
            initial_dict = {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}
        return self.process_coeff_dict(initial_dict)

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Same as ``cached_product``."""
        return self.cached_product(u, v, basis2)

    def __call__(self, x):
        """Build an element from a parabolic permutation/Lehmer list or an expression; raises ``ValueError``
        if the permutation is not parabolic for this ring's blocks.
        """
        genset = self.genset
        if not genset:
            genset = self.genset
        if isinstance(x, list) or isinstance(x, tuple):
            perm = Permutation(x)
            if not is_parabolic(perm, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {perm} is not")
            elem = self.from_dict({perm: self.domain.one})
        elif isinstance(x, Permutation):
            if not is_parabolic(x, self.parabolic_index):
                raise ValueError(f"Permutation must be parabolic: {x} is not")
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x
        else:
            elem = self.from_expr(x)
        return elem


def make_single_parabolic_quantum_basis(index_comp):
    """The `ParabolicQuantumSingleSchubertRing` in ``x`` for block sizes ``index_comp``."""
    return ParabolicQuantumSingleSchubertRing(GeneratingSet("x"), index_comp)


@cache
def QPSx(*args):
    """Cached `ParabolicQuantumSingleSchubertRing` for block sizes ``args``, e.g. ``QPSx(2, 1)([2, 1, 3])``."""
    return make_single_parabolic_quantum_basis(args)

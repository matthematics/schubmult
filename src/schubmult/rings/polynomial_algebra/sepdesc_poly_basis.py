from functools import cached_property

from schubmult.combinatorics.permutation import Permutation
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict

from ..schubert.schubert_ring import Sx
from .base_polynomial_basis import PolynomialBasis


class SepDescPolyBasis(PolynomialBasis):
    def __hash__(self):
        return hash((self.ring, self.k, "fatbacon"))

    @property
    def zero_monom(self):
        return self.as_key([[], []])

    def product(self, key1, key2, coeff=S.One):
        from .schubert_poly_basis import SchubertPolyBasis

        mnb = SchubertPolyBasis(ring=self.ring)
        left = self.transition(mnb)({key1: coeff})
        right = self.transition(mnb)({key2: S.One})

        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, mnb.transition(self)(mnb.product(key_schub_left, key_schub_right, v * v2)))
        return ret

    @property
    def genset(self):
        return self.ring.genset

    @property
    def k(self):
        return self._k

    def printing_term(self, k):
        return GenericPrintingTerm(k, "K")

    def is_key(self, x):
        return isinstance(x, list | tuple)

    def as_key(self, x):
        return (Permutation(x[0]), Permutation(x[1]), self.k)

    def __init__(self, k, ring=None):
        self._k = k
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    @cached_property
    def monomial_basis(self):
        from .monomial_basis import MonomialBasis

        return MonomialBasis(genset=self.ring.genset)

    def transition_schubert(self, dct, other_basis):
        out_ret = self.ring.zero
        for (k1, k2, _, __), v in dct.items():
            out_ret += self.ring.from_dict({k1: v}) * self.ring.from_dict({k2: S.One})
        return other_basis.attach_key(dict(out_ret))

    def transition_sepdesc(x, other_basis):  # noqa: ARG002
        return {x: S.One}

    def from_expr(self, expr):
        return self.attach_key(self.ring.from_expr(expr))

    def transition(self, other_basis):
        from .elem_sym_poly_basis import ElemSymPolyBasis
        from .monomial_basis import MonomialBasis
        from .schubert_poly_basis import SchubertPolyBasis

        if isinstance(other_basis, SepDescPolyBasis):
            return lambda x: self.transition_sepdesc(x, other_basis)
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: self.transition_schubert(x, other_basis)
        if isinstance(other_basis, MonomialBasis):
            bonky_basis = SchubertPolyBasis(self.ring)
            return lambda x: bonky_basis.transition(other_basis)(self.transition_schubert(x, bonky_basis))
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: other_basis.transition_elementary(self.transition_schubert(x, other_basis))
        return None

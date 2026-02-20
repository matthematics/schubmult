from functools import cached_property

from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S

from ..schubert.schubert_ring import Sx
from .base_polynomial_basis import PolynomialBasis


class ElemSymPolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return ((*x,),)

    def printing_term(self, k):
        return GenericPrintingTerm(k, "E")

    def __hash__(self):
        return hash((self.ring, "dangit_bobjo"))

    @cached_property
    def monomial_basis(self):
        from .monomial_basis import MonomialBasis

        return MonomialBasis(genset=self.ring.genset)

    @property
    def genset(self):
        return self.ring.genset

    def __init__(self, ring=None):
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    def transition_schubert(self, dct):
        from schubmult.abc import e

        res = self.ring.zero
        for (k, n), v in dct.items():
            to_add = self.ring.one
            for i, a in enumerate(k[:n]):
                to_add = self.ring.elem_mul(to_add, e(a, i + 1, self.ring.genset[1:]))
            for a in enumerate(k[n:]):
                to_add = self.ring.elem_mul(to_add, e(a, n, self.ring.genset[1:]))
            res += v * to_add
        return res

    def transition_monomial(self, dct):
        from schubmult.abc import e
        from schubmult.symbolic import expand_func

        res = S.Zero
        for (k, n), v in dct.items():
            to_add = S.One
            for i, a in enumerate(k[:n]):
                to_add *= expand_func(e(a, i + 1, self.ring.genset[1:]))
            for a in k[n:]:
                to_add *= expand_func(e(a, n, self.ring.genset[1:]))
            res += v * to_add
        return res

    def transition(self, other_basis):
        from .monomial_basis import MonomialBasis
        from .schubert_poly_basis import SchubertPolyBasis

        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: other_basis.attach_key(self.transition_schubert(x))
        if isinstance(other_basis, MonomialBasis):
            from schubmult.symbolic.poly.variables import genset_dict_from_expr

            return lambda x: other_basis.attach_key(genset_dict_from_expr(self.transition_monomial(x), other_basis.genset))
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: x
        return None

    @property
    def zero_monom(self):
        return self.as_key([])

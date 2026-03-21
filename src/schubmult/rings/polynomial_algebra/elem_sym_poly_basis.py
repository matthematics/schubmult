from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S

from ..schubert.schubert_ring import SingleSchubertRing
from .base_polynomial_basis import PolynomialBasis


class ElemSymPolyBasis(PolynomialBasis):
    """Elementary symmetric polynomial basis.

    Keys are tuples encoding products of elementary symmetric polynomials
    e_k(x_1, ..., x_n). Each key specifies degrees and variable counts
    for the elementary symmetric factors.
    """
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return ((*x,),)

    def printing_term(self, k):
        return GenericPrintingTerm(k, "E")

    def __hash__(self):
        return hash((self.ring, "dangit_bobjo"))

    def __init__(self, genset):
        from .monomial_basis import MonomialBasis
        self.ring = SingleSchubertRing(genset)
        super().__init__(genset=genset)
        self._monomial_basis = MonomialBasis(genset=genset)

    def transition_schubert(self, dct):
        """Transition from elementary symmetric basis to Schubert basis."""
        from schubmult.abc import e

        res = {}
        for (k, n), v in dct.items():
            to_add = self.ring.one
            for i, a in enumerate(k[:n]):
                if a == 0:
                    continue
                to_add = self.ring.elem_mul(to_add, e(a, i + 1, self.ring.genset[1:]))
            for i, a in enumerate(k[n:]):
                if a == 0:
                    continue
                to_add = self.ring.elem_mul(to_add, e(a, n, self.ring.genset[1:]))
            for perm, coeff in to_add.items():
                res[(perm, n)] = res.get((perm, n), S.Zero) + v * coeff
        return res

    def transition_monomial(self, dct):
        """Transition from elementary symmetric basis to monomial basis."""
        from schubmult.abc import e
        from schubmult.symbolic import expand_func
        from schubmult.utils.perm_utils import add_perm_dict_with_coeff

        res = {}
        for (k, n), v in dct.items():
            to_add = S.One
            for i, a in enumerate(k[:n]):
                to_add *= expand_func(e(a, i + 1, self.ring.genset[1:]))
            for i, a in enumerate(k[n:]):
                to_add *= expand_func(e(a, n, self.ring.genset[1:]))
            res = add_perm_dict_with_coeff(res, self.monomial_basis.from_expr(to_add, length=n), coeff=v)
        return res

    def transition(self, other_basis):
        """Return a transition function from this basis to *other_basis*."""
        from .monomial_basis import MonomialBasis
        from .schubert_poly_basis import SchubertPolyBasis

        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: other_basis.attach_key(self.transition_schubert(x))
        if isinstance(other_basis, MonomialBasis):
            return lambda x: other_basis.attach_key(self.transition_monomial(x))
        if isinstance(other_basis, ElemSymPolyBasis):
            return lambda x: x
        return None

    @property
    def zero_monom(self):
        return self.as_key([])

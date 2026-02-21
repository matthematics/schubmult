from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import Add, Mul, S
from schubmult.utils.perm_utils import add_perm_dict

from .base_polynomial_basis import PolynomialBasis


class MonomialBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(str(self.expand_monom(k)), "")

    def coproduct(self, key):
        result_dict = {}
        key = self.as_key(key)
        for i in range(len(key) + 1):
            result_dict[(key[:i], key[i:])] = S.One
        return result_dict

    @property
    def monomial_basis(self):
        return self

    @property
    def genset(self):
        return self._genset

    def __init__(self, genset):
        self._genset = genset

    def product(self, key1, key2, coeff=S.One):
        if len(key1) != len(key2):
            return {}
        return {tuple(a + b for a, b in zip(key1, key2)): coeff}

    def expand_monom(self, monom):
        return Mul(*[self.genset[i + 1] ** monom[i] for i in range(len(monom))])

    def expand(self, dct):
        return Add(*[v * self.expand_monom(k) for k, v in dct.items()])

    def transition_slide(self, other_basis, dct):
        ret = {}
        for k, v in dct.items():
            ret = add_perm_dict(ret, self.transition_slide_monom(other_basis, k, coeff=v))
        return ret

    def transition_slide_monom(self, other_basis, monom, coeff=S.One):
        from .slide_poly_basis import SlidePolyBasis

        def domkey(comp):
            return tuple([sum(comp[:i]) for i in range(1, len(comp))])

        if not isinstance(other_basis, SlidePolyBasis):
            return None

        res = {monom: coeff}
        ret = {}
        while len(res) > 0:
            key = next(iter(sorted(res.keys(), key=lambda x: domkey(x))))
            val = res[key]
            dct2 = {k: v * val for k, v in other_basis.to_monoms(key).items()}
            ret[key] = ret.get(key, 0) + val
            for k, v in dct2.items():
                res[k] = res.get(k, 0) - v
                if res[k] == 0:
                    del res[k]
        return ret

    def transition(self, other_basis):
        from .elem_sym_poly_basis import ElemSymPolyBasis
        from .schubert_poly_basis import SchubertPolyBasis
        from .sepdesc_poly_basis import SepDescPolyBasis
        from .slide_poly_basis import SlidePolyBasis

        if isinstance(other_basis, MonomialBasis):
            return lambda x: other_basis.attach_key(x)
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: other_basis.attach_key(other_basis.ring.from_expr(Add(*[v * self.expand_monom(k) for k, v in x.items()])))
        if isinstance(other_basis, SlidePolyBasis):
            return lambda x: self.transition_slide(other_basis, x)
        if isinstance(other_basis, SepDescPolyBasis):
            bonky_basis = SchubertPolyBasis(ring=other_basis.ring)
            return lambda x: other_basis.attach_key(bonky_basis.transition(other_basis)(bonky_basis.attach_key(bonky_basis.ring.from_expr(Add(*[v * self.expand_monom(k) for k, v in x.items()])))))
        if isinstance(other_basis, ElemSymPolyBasis):
            spb = SchubertPolyBasis(ring=other_basis.ring)
            return lambda x: spb.transition(other_basis)(self.transition(spb)(x))
        return None

    def from_expr(self, expr):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        return genset_dict_from_expr(expr, self.genset)

    @property
    def zero_monom(self):
        return self.as_key([])

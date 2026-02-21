from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import Add, S, expand_seq
from schubmult.utils.perm_utils import add_perm_dict_with_coeff

from .base_polynomial_basis import PolynomialBasis


def _slide_polynomial(comp, genset):
    if len(comp) == 0 or all(c == 0 for c in comp):
        return S.One

    if all(c > 0 for c in comp):
        return expand_seq(comp, genset)

    ret = S.Zero

    min_nonzero = min(i for i, c in enumerate(comp) if c > 0)

    power = comp[min_nonzero]

    for i in range(min_nonzero + 1):
        ret += genset[min_nonzero + 1 - i] ** power * _slide_polynomial((0,) * i + comp[min_nonzero + 1 :], genset[min_nonzero + 1 - i :])

    return ret


class MonomialSlidePolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"MSlide{k}", "")

    # def coproduct(self, key):
    #     result_dict = {}
    #     key = self.as_key(key)
    #     for i in range(len(key) + 1):
    #         result_dict[(key[:i], key[i:])] = S.One
    #     return result_dict

    @property
    def monomial_basis(self):
        return self._monomial_basis

    @property
    def genset(self):
        return self._genset

    def __init__(self, genset):
        from .monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        def pad_tuple(tup, length):
            return (*tup, *(0,) * (length - len(tup)))

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_slide_polynomial(key, self.genset), self.genset).items()}
        return dct

    def transition_monomial(self, dct):
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def expand(self, dct):
        return Add(*[v * self.to_monoms(k) for k, v in dct.items()])

    def transition(self, other_basis):
        if isinstance(other_basis, self.monomial_basis.__class__):
            return self.transition_monomial
        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def from_expr(self, expr):
        dct = self.monomial_basis.from_expr(expr)
        return self.monomial_basis.transition_slide(self, dct)

    @property
    def zero_monom(self):
        return self.as_key([])

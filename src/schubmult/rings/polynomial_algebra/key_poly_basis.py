from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.plactic import Plactic
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import add_perm_dict_with_coeff
from schubmult.utils.tuple_utils import pad_tuple

from ..printing import GenericPrintingTerm
from .base_polynomial_basis import PolynomialBasis

"""
Fundamental slide polynomial basis for Schubert calculus.

This module implements the fundamental slide polynomial basis, which provides
an alternative basis for expressing Schubert polynomials and their products.
"""

def traverse_demaz(pl, w):
    stack = [(pl, 0)]
    word = list(reversed(w.code_word))
    while stack:
        current, index = stack.pop()
        yield current
        for index2 in range(index, len(word)):
            desc2 = word[index2]
            next_elem = current.lowering_operator(desc2)
            if next_elem is not None:
                stack.append((next_elem, index2))




def _key_polynomial(comp, genset):
    w = Permutation.sorting_perm(comp, reverse=True)
    hw = sorted(comp, reverse=True)
    top_pl = Plactic.yamanouchi(hw)
    visited = set()
    ret = S.Zero
    for pl in traverse_demaz(top_pl, w):
        if pl in visited:
            continue
        visited.add(pl)
        ret += expand_seq(pl.crystal_weight, genset)
    return ret

class KeyPolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"KeyPoly{k}", "")

    @classmethod
    def dual_basis(cls):
        from ..free_algebra.key_basis import KeyBasis
        return KeyBasis

    def __init__(self, genset):
        from .monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_key_polynomial(key, self.genset), self.genset).items()}
        return dct

    def expand(self, dct):
        return sum([v * _key_polynomial(key, self.genset) for key, v in dct.items()])

    def transition_monomial(self, dct):
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition(self, other_basis):
        from .monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def from_expr(self, expr, length=None):
        _ = length
        # dct = self.monomial_basis.from_expr(expr)
        # return self.monomial_basis.transition_slide(dct, self)
        try:
            from schubmult.symbolic import sympify
            return {self.zero_monom: sympify(expr)}
        except Exception:
            pass
        raise NotImplementedError("Direct expression parsing not implemented for FundamentalSlidePolyBasis yet")

    @property
    def zero_monom(self):
        return self.as_key([])

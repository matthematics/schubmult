from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor
from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict_with_coeff


def _forest_polynomial_from_indfor(indfor, genset):
    if len(indfor) == 0:
        return S.One

    ret = S.One

    for root in indfor:
        ret *= _forest_polynomial_from_root(root, genset)
    return ret


def _forest_polynomial_from_root(root, genset, min_value=1):
    if root.rho == 0:
        return S.One
    if root.left is None and root.right is None:
        return sum([genset[i] for i in range(min_value, root.rho + 1)])
    ret = S.Zero
    for val in range(min_value, root.rho + 1):
        term = genset[val]
        if root.left is not None:
            term *= _forest_polynomial_from_root(root.left, genset, val)
        if root.right is not None:
            term *= _forest_polynomial_from_root(root.right, genset, val + 1)
        ret += term
    return ret


class ForestPolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"For{k}", "")

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
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        def pad_tuple(tup, length):
            return (*tup, *(0,) * (length - len(tup)))

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_forest_polynomial_from_indfor(weak_composition_to_indfor(key), self.genset), self.genset).items()}
        return dct

    def expand(self, dct):
        return sum([v * _forest_polynomial_from_indfor(weak_composition_to_indfor(key), self.genset) for key, v in dct.items()])

    def transition_monomial(self, dct):
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition(self, other_basis):
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def from_expr(self, expr):
        dct = self.monomial_basis.from_expr(expr)
        return self.monomial_basis.transition_slide(self, dct)

    @property
    def zero_monom(self):
        return self.as_key([])


if __name__ == "__main__":
    from symengine import expand

    from schubmult import PolynomialAlgebra, Sx

    For = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    print(expand(For([0, 2, 0, 1]).expand()))

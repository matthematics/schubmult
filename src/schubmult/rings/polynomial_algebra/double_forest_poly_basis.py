from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import expand
from schubmult.utils.perm_utils import add_perm_dict_with_coeff
from schubmult.utils.tuple_utils import pad_tuple


class DoubleForestPolyBasis(PolynomialBasis):
    """Abstract double forest polynomial basis.

    Keys are weak compositions indexing double forest basis elements DF[key],
    with polynomial coefficients in a second generating set (equivariant vars).
    """

    def __init__(self, genset, coeff_genset, n=None):
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        super().__init__(genset=genset)
        self.coeff_genset = coeff_genset
        self._n = n
        self._basis_poly_cache = {}
        self._basis_forest_cache = {}
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    @property
    def zero_monom(self):
        return self.as_key([])

    def printing_term(self, k):
        return GenericPrintingTerm(f"DForest{k}", "")

    def _x_gen(self, i):
        return self.genset[i]

    def _t_gen(self, i):
        return self.coeff_genset[i]

    def basis_polynomial(self, key):
        """Expand one double-forest basis key as a polynomial in x,t."""
        from schubmult.combinatorics.double_forest import double_forest_polynomial

        k = self.as_key(key)
        if k not in self._basis_poly_cache:
            self._basis_poly_cache[k] = expand(double_forest_polynomial(k, self.genset, self.coeff_genset, n=self._n))
        return self._basis_poly_cache[k]

    def basis_forest_expansion(self, key, length):
        """Expand one double-forest basis key into ForestPolyBasis in x."""
        from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra

        k = self.as_key(key)
        cache_key = (k, length)
        if cache_key not in self._basis_forest_cache:
            forest_ring = PolynomialAlgebra(ForestPolyBasis(self.genset))
            self._basis_forest_cache[cache_key] = {
                tuple(kk): expand(vv) for kk, vv in forest_ring.from_expr(self.basis_polynomial(k), length=length).items() if vv != 0
            }
        return self._basis_forest_cache[cache_key]

    def _clean_dict(self, dct):
        return {tuple(k): expand(v) for k, v in dct.items() if expand(v) != 0}

    def _pad_and_combine(self, dct, length):
        out = {}
        for k, v in dct.items():
            kk = pad_tuple(tuple(k), length)
            out[kk] = expand(out.get(kk, 0) + v)
        return self._clean_dict(out)

    def _straighten_from_forest_dict(self, forest_dct, length):
        """Peel top x-degree terms in Forest basis until residual vanishes."""
        residual = self._pad_and_combine(forest_dct, length)
        out = {}

        for _ in range(10000):
            if not residual:
                break
            top_degree = max(sum(k) for k in residual)
            top_keys = [k for k in residual if sum(k) == top_degree]

            progressed = False
            for key in sorted(top_keys, reverse=True):
                coeff = expand(residual.get(key, 0))
                if coeff == 0:
                    continue
                out[key] = expand(out.get(key, 0) + coeff)

                df_as_forest = self.basis_forest_expansion(key, length)
                residual = add_perm_dict_with_coeff(residual, df_as_forest, coeff=-coeff)
                residual = self._pad_and_combine(residual, length)
                progressed = True

            if not progressed:
                raise RuntimeError("Double forest straightening stalled")

        if residual:
            raise RuntimeError(f"Double forest straightening did not terminate: {residual}")

        return self._clean_dict(out)

    def from_forest_dict(self, forest_dct, length):
        return self._straighten_from_forest_dict(forest_dct, length=length)

    def expand(self, dct):
        return sum([v * self.basis_polynomial(k) for k, v in dct.items()])

    def transition_monomial(self, dct):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        res = {}
        for k, v in dct.items():
            poly = self.basis_polynomial(k)
            monoms = genset_dict_from_expr(poly, self.genset, length=len(k))
            res = add_perm_dict_with_coeff(res, monoms, coeff=v)
        return res

    def transition_forest(self, dct):
        res = {}
        for k, v in dct.items():
            len_k = len(k)
            as_forest = self.basis_forest_expansion(k, length=len_k)
            res = add_perm_dict_with_coeff(res, as_forest, coeff=v)
        return res

    def transition(self, other_basis):
        from schubmult.rings.polynomial_algebra.forest_poly_basis import ForestPolyBasis
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial
        if isinstance(other_basis, ForestPolyBasis):
            return self.transition_forest

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

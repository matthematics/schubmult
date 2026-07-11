from schubmult.combinatorics.indexed_forests import grove_polynomial
from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import Symbol
from schubmult.utils.perm_utils import add_perm_dict_with_coeff

"""
Grove polynomial basis for Schubert calculus.

Grove polynomials are the set-valued (``beta``-deformed) analog of forest
polynomials. Keys are weak compositions encoding indexed forests; the grove
polynomial sums over set-valued labelings of the forest, weighting each extra
label by ``beta`` (see :func:`grove_polynomial`).
"""


class GrovePolyBasis(PolynomialBasis):
    """Grove polynomial basis.

    Keys are weak compositions encoding indexed forests. Grove polynomials are
    the ``beta``-deformed (set-valued) forest polynomials, computed by summing
    over set-valued labelings of the corresponding forest structure.
    """

    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"Grove{k}", "")

    def __init__(self, genset, beta=Symbol("beta")):
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        super().__init__(genset=genset)
        self._beta = beta
        self._monomial_basis = MonomialBasis(genset=self.genset)

    @property
    def beta(self):
        return self._beta

    def _grove_poly(self, key):
        return grove_polynomial(key, self.genset, beta=self.beta)

    def to_monoms(self, key):
        """Expand a grove key into a dict of monomial exponent tuples."""
        from schubmult import WCGraph

        dct = {}
        for wc in WCGraph.grove_wcs(key):
            dct[wc.length_vector] = dct.get(wc.length_vector, 0) + 1
        #dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(self._grove_poly(key), self.genset).items()}
        return dct

    def expand(self, dct):
        """Expand a grove basis dict into a symbolic polynomial expression."""
        return sum([v * self._grove_poly(key) for key, v in dct.items()])

    def transition_monomial(self, dct):
        """Transition from grove basis to monomial basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition(self, other_basis):
        """Return a transition function from grove basis to *other_basis*."""
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, GrovePolyBasis):
            return lambda x: x
        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    @property
    def zero_monom(self):
        return self.as_key([])


if __name__ == "__main__":
    from symengine import expand

    from schubmult.abc import x

    basis = GrovePolyBasis(x)
    print(expand(basis.expand({(0, 2, 3): 1})))

from schubmult.combinatorics.permutation import uncode
from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.utils.perm_utils import add_perm_dict_with_coeff

"""
Fundamental slide polynomial basis for Schubert calculus.

This module implements the fundamental slide polynomial basis, which provides
an alternative basis for expressing Schubert polynomials and their products.
"""


# def _glide_polynomial(key, genset):
#     """Compute the glide polynomial for a weak composition key."""
#     from schubmult import WCGraph

#     dct = {}
#     for wc in WCGraph.glide_wcs(key):
#         dct[wc.length_vector] = dct.get(wc.length_vector, 0) + 1
#     return PA.from_dict(dct, genset=genset)


class GlidePolyBasis(PolynomialBasis):
    """Glide polynomial basis.

    Keys are weak compositions. Glide polynomials provide a
    basis that refines Grothendieck polynomials and coarsens monomials, with
    an efficient combinatorial product rule.
    """

    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"GlidePoly{k}", "")

    def __init__(self, genset):
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=genset)

    def to_monoms(self, key):
        """Expand a glide key into a dict of monomial exponent tuples."""
        from schubmult.combinatorics.wc_graph import WCGraph
        key = tuple(key)
        the_set = {wc for wc in WCGraph.all_wc_graphs(uncode(key), len(key)) if wc.dst == WCGraph.principal_wc(uncode(key), len(key))}
        dct = {}
        for wc in the_set:
            dct[wc.length_vector] = dct.get(wc.length_vector, 0) + 1
        return dct

    @classmethod
    def dual_basis(cls):
        """Return the dual free algebra basis class (:class:`GlideBasis`)."""
        from ..free_algebra.fundamental_slide_basis import GlideBasis

        return GlideBasis

    def expand(self, dct):
        """Expand a glide basis dict into a symbolic polynomial expression."""
        return sum([v * self._monomial_basis.from_dict(self.to_monoms(key)).expand() for key, v in dct.items()])

    def transition_monomial(self, dct):
        """Transition from glide basis to monomial basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition(self, other_basis):
        """Return a transition function from glide basis to *other_basis*."""
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self._monomial_basis.transition(other_basis), self.transition_monomial(x))

    # def product(self, key1, key2, coeff=S.One):
    #     """Multiply two glide keys using the glide product rule."""
    #     return {c: v * coeff for c, v in glide_product(key1, key2).items()}

    @property
    def zero_monom(self):
        return self.as_key([])


if __name__ == "__main__":
    from schubmult import PolynomialAlgebra, Sx

    FSlide = PolynomialAlgebra(GlidePolyBasis(Sx.genset))
    print(FSlide(0, 2, 0, 3) * FSlide(1, 0, 0, 1))

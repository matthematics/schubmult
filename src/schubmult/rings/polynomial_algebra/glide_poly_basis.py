from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict_with_coeff

"""
Fundamental slide polynomial basis for Schubert calculus.

This module implements the fundamental slide polynomial basis, which provides
an alternative basis for expressing Schubert polynomials and their products.
"""


def glide_monomials(key):
    """Monomial expansion of the glide polynomial :math:`\\mathcal{G}_{key}`.

    Returns a dict mapping each exponent tuple ``v`` (a weak composition of the
    same length as ``key``) to the integer coefficient of the monomial
    :math:`x^v`, i.e. the number of glides of ``key`` with weight ``v``. The
    corresponding power of ``beta`` for the weight ``v`` is
    ``sum(v) - sum(key)`` (the excess), which is constant across all glides of a
    given weight, so it need not be stored explicitly.
    """
    from schubmult.combinatorics.wc_graph import WCGraph

    key = tuple(key)
    perm = uncode(key)
    principal = WCGraph.principal_wc(perm, len(key))
    dct = {}
    for wc in WCGraph.all_wc_graphs(perm, len(key)):
        if wc.dst == principal:
            v = wc.length_vector
            dct[v] = dct.get(v, 0) + 1
    return dct


def glide_product(key1, key2):
    """Structure constants for a product of two glide polynomials.

    Implements the Littlewood-Richardson rule of O. Pechenik and D. Searles,
    "Decompositions of Grothendieck Polynomials" (arXiv:1611.02545), Theorem
    4.9, which expands the product of the glide polynomials indexed by the weak
    compositions ``key1`` and ``key2`` in the glide basis:

    .. math::

        \\mathcal{G}_a \\, \\mathcal{G}_b
            = \\sum_c \\beta^{|c| - |a| - |b|} \\, g_{a,b}^{c} \\, \\mathcal{G}_c .

    Rather than enumerating the genomic shuffle set directly, we compute the
    (uniquely determined) coefficients by expanding the product in monomials and
    straightening into the glide basis with the leading-term algorithm from the
    proof that the glide polynomials form a basis (Theorem 2.6). Because the
    excess of a glide of ``v`` equals ``sum(v) - sum(index)``, the power of
    ``beta`` is recovered from the total degree and only the positive integer
    multiplicities :math:`g_{a,b}^{c}` are returned.

    Both compositions are padded with trailing zeros to a common length ``n``;
    every key ``c`` in the returned dict is a weak composition of length ``n``.
    """
    a = tuple(key1)
    b = tuple(key2)
    n = max(len(a), len(b))
    a = a + (0,) * (n - len(a))
    b = b + (0,) * (n - len(b))

    # Monomial expansion of the product G_a * G_b as {weight: coefficient}.
    ma = glide_monomials(a)
    mb = glide_monomials(b)
    product = {}
    for va, ca in ma.items():
        for vb, cb in mb.items():
            v = tuple(x + y for x, y in zip(va, vb))
            product[v] = product.get(v, 0) + ca * cb

    def leading_key(v):
        # Total order of arXiv:1611.02545, Thm 2.6: strictly more zeros is
        # larger; ties are broken by reverse-lexicographic order (ordinary
        # lexicographic order on the reversed vector).
        return (v.count(0), tuple(reversed(v)))

    result = {}
    while product:
        v_star = max(product, key=leading_key)
        # The leading monomial of G_{v_star} is x^{v_star} with coefficient 1,
        # so its coefficient in the product is the structure constant.
        coeff = product.pop(v_star)
        if coeff == 0:
            continue
        result[v_star] = coeff
        for w, cnt in glide_monomials(v_star).items():
            if w == v_star:
                continue
            new = product.get(w, 0) - coeff * cnt
            if new == 0:
                product.pop(w, None)
            else:
                product[w] = new
    return result


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
        return glide_monomials(key)

    @classmethod
    def dual_basis(cls):
        """Return the dual free algebra basis class (:class:`GlideBasis`)."""
        from ..free_algebra.fundamental_slide_basis import GlideBasis

        return GlideBasis

    def expand(self, dct):
        """Expand a glide basis dict into a symbolic polynomial expression."""
        return sum([v * self._monomial_basis.expand(self.to_monoms(key)) for key, v in dct.items()])

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

    @cache
    def product(self, key1, key2, coeff=S.One):
        """Multiply two glide keys using the glide product rule."""
        return {c: v * coeff for c, v in glide_product(key1, key2).items()}

    @property
    def zero_monom(self):
        return self.as_key([])


if __name__ == "__main__":
    from schubmult import PolynomialAlgebra, Sx

    FSlide = PolynomialAlgebra(GlidePolyBasis(Sx.genset))
    print(FSlide(0, 2, 0, 3) * FSlide(1, 0, 0, 1))

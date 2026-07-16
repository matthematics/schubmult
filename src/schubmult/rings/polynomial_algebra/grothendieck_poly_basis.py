from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict, add_perm_dict_with_coeff
from schubmult.utils.tuple_utils import pad_tuple

from .base_polynomial_basis import PolynomialBasis


class GrothendieckPolyBasis(PolynomialBasis):
    """Grothendieck polynomial basis.

    Keys are ``(Permutation, length)`` pairs. Grothendieck polynomials form
    the canonical basis for the polynomial algebra in Grothendieck calculus,
    dual to the :class:`GrothendieckBasis` of the free algebra.
    """
    def __hash__(self):
        return hash(("schoobooponk", self.ring))

    # def coproduct(self, key):
    #     """Compute the coproduct of a Grothendieck key by splitting variable sets."""
    #     length = key[1]
    #     res = {}
    #     for len0 in range(length + 1):
    #         cprd = dict(self.ring(key[0]).coproduct(*list(range(1, len0 + 1))))
    #         res = add_perm_dict(res, {((k1, len0), (k2, length - len0)): v for (k1, k2), v in cprd.items()})
    #     return res

    def printing_term(self, k):
        return self.ring.printing_term(k)

    def is_key(self, x):
        return isinstance(x, Permutation) or (isinstance(x, list | tuple) and len(x) == 2 and isinstance(x[0], Permutation) and isinstance(x[1], int))

    def as_key(self, x):
        if isinstance(x, Permutation):
            return (x, len(x.trimcode))
        return tuple(x)

    def __init__(self, genset=None, ring=None):
        from .monomial_basis import MonomialBasis

        if ring is not None:
            self.ring = ring
            working_genset = self.ring.genset
        elif hasattr(genset, "from_dict") and hasattr(genset, "genset"):
            self.ring = genset
            working_genset = self.ring.genset
        else:
            if genset is None:
                from ..schubert.schubert_ring import Sx

                self.ring = Sx
                working_genset = self.ring.genset
            else:
                from ..schubert.grothendieck_ring import GrothendieckRing
                self.ring = GrothendieckRing(genset, beta=1)
                working_genset = genset

        super().__init__(genset=working_genset)
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def product(self, key1, key2, coeff=S.One):
        """Multiply two Grothendieck keys using the Grothendieck ring multiplication."""
        if key1[1] != key2[1]:
            return {}

        def pair_length(dct, length):
            return {(k, length): v for k, v in dct.items()}

        return pair_length(self.ring.mul(self.ring.from_dict({key1[0]: coeff}), self.ring(key2[0])), key1[1])

    def transition_schubert(self, dct):
        """Transition from Grothendieck basis to separated descents basis."""
        from schubmult import WCGraph
        ret = {}
        for (perm, length), coeff in dct.items():
            for perm2, coeff2 in WCGraph.groth_to_schub(perm, self.ring._beta).items():
                #yield (perm2, length), coeff * coeff2
                ret[(perm2, length)] = ret.get((perm2, length), S.Zero) + coeff * coeff2
        return ret


    # def transition_slide(self, dct, other_basis):
    #     from schubmult.abc import x
    #     from schubmult.symbolic import expand_func

    #     if not isinstance(other_basis, GrothendieckPolyBasis):

    @property
    def zero_monom(self):
        return (Permutation([]), 0)

    def transition_glide_key(self, key):
        """Decompose a Grothendieck polynomial into glide polynomials via omega insertion on RC-graphs."""
        from schubmult import WCGraph

        dct = {}
        for wc in WCGraph.all_wc_graphs(key[0], key[1]):
            if wc.is_quasi_yamanouchi:
                dct[wc.length_vector] = dct.get(wc.length_vector, S.Zero) + S.One
        return dct

    def transition_glide(self, dct):
        """Transition a Grothendieck dict to the glide polynomial basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.transition_glide_key(k), coeff=v)
        return res

    def to_monoms(self, key):
        """Expand a Grothendieck key into a dict of monomial exponent tuples."""
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        return {pad_tuple(k, key[1]): v for k, v in genset_dict_from_expr(self.ring.from_dict({key[0]: S.One}).as_polynomial(), self.genset).items()}

    @classmethod
    def dual_basis(cls):
        """Return the dual free algebra basis class (:class:`GrothendieckBasis`)."""
        from ..free_algebra.schubert_basis import GrothendieckBasis
        return GrothendieckBasis

    def transition_grove_key(self, key):
        """Decompose a Grothendieck polynomial into grove polynomials via omega insertion on RC-graphs."""
        from schubmult import WCGraph

        dct = {}
        for wc in WCGraph.all_wc_graphs(key[0], key[1]):
            if wc.grove_weight == wc.length_vector:
                keykey = wc.grove_weight
                dct[keykey] = dct.get(keykey, S.Zero) + S.One
        return dct

    def transition_grove(self, dct):
        """Transition a Grothendieck dict to the grove polynomial basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.transition_grove_key(k), coeff=v)
        return res

    def transition(self, other_basis):
        """Return a transition function from Grothendieck basis to *other_basis*."""
        from .glide_poly_basis import GlidePolyBasis
        from .grove_poly_basis import GrovePolyBasis
        from .monomial_basis import MonomialBasis
        from .schubert_poly_basis import SchubertPolyBasis


        if isinstance(other_basis, GrothendieckPolyBasis):
            return lambda x: dict(x)
        if isinstance(other_basis, GlidePolyBasis):
            return lambda x: self.transition_glide(x)
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: self.transition_schubert(x)
        if isinstance(other_basis, MonomialBasis):
            def sum_dct(*dcts):
                res = {}
                for dct in dcts:
                    res = add_perm_dict(res, dct)
                return res

            def mul_dict(dct, coeff):
                return {k: v * coeff for k, v in dct.items()}

            return lambda x: sum_dct(*[mul_dict(self.to_monoms(k), v) for k, v in x.items()])
        if isinstance(other_basis, GrovePolyBasis):
            return lambda x: self.transition_grove(x)

        return lambda x: PolynomialBasis.compose_transition(SchubertPolyBasis(self.monomial_basis.genset).transition(other_basis), self.transition_schubert(x))

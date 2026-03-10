from functools import cache

from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.combinatorial.schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement
from schubmult.symbolic import S

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class DualRCGraphRingElement(SchubertMonomialRingElement):
    """
    DualRCGraphRing elements are linear combinations of RCGraph basis elements.

    The product % is the polynomial product. Currently only defined when the right side
    is a dominant RC graph.

    The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
    to do this directly.

    The product * is well defined for any pair of RC graphs and is the dual product.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------

    def divdiff_perm(self, perm):
        """
        Apply divided difference operator for `perm` to self.
        Linear extension of RCGraph.divdiff_perm.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc_set = rc_graph.divdiff_perm(perm)
            for new_rc in new_rc_set:
                res += coeff * self.ring(new_rc)
        return res

    def divdiff(self, *seq):
        """
        Sequential divided difference operators.
        """
        total = self
        for i in reversed(seq):
            res = self.ring.zero
            for rc_graph, coeff in total.items():
                new_rc_set = rc_graph.divdiff_desc(i)
                for new_rc in new_rc_set:
                    res += coeff * self.ring(new_rc.resize(len(rc_graph)))
            total = res
        return total

    # def coproduct(self):
    #     """
    #     Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
    #     and induces the mul product.
    #     """
    #     tring = self.ring @ self.ring
    #     res = tring.zero
    #     for rc_graph, coeff in self.items():
    #         for i in range(len(rc_graph) + 1):
    #             rc1, rc2 = rc_graph.vertical_cut(i)
    #             res += tring.from_dict({(rc1, rc2): coeff})
    #     return res

    def almosteq(self, other):
        return all(v == 0 for v in (self - other).values())

    def shiftup(self, k):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.shiftup(k))
        return res

    def prepend(self, k):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.prepend(k))
        return res

    def zero_out_last_row(self):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.zero_out_last_row())
        return res

    def resize(self, n):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.resize(n))
        return res

    def clip(self, n):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.vertical_cut(n)[0])
        return res

    def transpose(self, length):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.transpose(length))
        return res


class DualRCGraphRing(SchubertMonomialRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = DualRCGraphRing._id
        DualRCGraphRing._id += 1
        self.dtype = type("DualRCGraphRingElement", (DualRCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtwtrwystoa", self._ID))

    def elem_sym(self, descent, weight):
        from schubmult import uncode

        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(uncode([0] * (descent - sum(weight)) + [1] * sum(weight)), len(weight), weight=weight), 1))

    def full_elem_sym(self, degree, descent, length):
        from schubmult import uncode

        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(uncode([0] * (descent - degree) + [1] * degree), length), 1))

    @property
    def zero_monom(self):
        return RCGraph([])

    # def from_free_algebra_element(self, elem):
    #     wordelem = elem.change_basis(WordBasis)
    #     result = self.zero
    #     for word, coeff in wordelem.items():
    #         res = self(RCGraph([]))
    #         for a in reversed(word):
    #             res = self(RCGraph.one_row(a)) * res
    #         result += coeff * res
    #     return result

    def schub(self, perm, n=None):
        """
        Return the DualRCGraphRing element corresponding to the Schubert polynomial
        indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).
        """
        if n is None:
            n = len(perm.trimcode)
        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, n), S.One))

    # def weight_coproduct(self, elem):

    @cache
    def coproduct_on_basis(self, elem):
        tring = self.ring @ self.ring
        res = tring.zero
        for i in range(len(elem) + 1):
            rc1, rc2 = elem.vertical_cut(i)
            res += tring.from_dict({(rc1, rc2): S.One})
        return res

    @property
    def one(self):
        # Define the "one" element for DualRCGraphRing
        identity_graph = RCGraph()
        return self.from_dict({identity_graph: 1})

    def mul(self, elem1, elem2):
        from .rc_graph_ring import RCGraphRing
        r = RCGraphRing()
        # Define the multiplication for DualRCGraphRing

        if isinstance(elem2, DualRCGraphRingElement):
            result = self.zero
            for a, coeff_a in elem1.items():
                for b, coeff_b in elem2.items():
                    if len(a) != len(b):
                        continue
                    result += coeff_a * coeff_b * self.from_dict(r(a)%r(b))

            return result
        try:
            result_dict = {k: v * elem2 for k, v in elem1.items()}
            return self.from_dict(result_dict)
        except Exception:
            raise NotImplementedError("Multiplication with non-DualRCGraphRingElement not implemented")

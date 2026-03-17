from functools import cache

from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRingElement
from schubmult.rings.combinatorial.schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement
from schubmult.symbolic import S

from ..polynomial_algebra import PA, SchubertPolyBasis

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

    def __rmul__(self, other):
        if isinstance(other, RCGraphRingElement):
            return self.rc_ring_left_act(other, self)
        return super().__rmul__(other)

    def __mul__(self, other):
        if isinstance(other, RCGraphRingElement):
            return self.rc_ring_right_act(self, other)
        return super().__mul__(other)

    def polynomial_coaction(self, left=True):
        res = PA.zero @ self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring.polynomial_coaction(rc_graph, left=left)
        return res


class DualRCGraphRing(SchubertMonomialRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = DualRCGraphRing._id
        DualRCGraphRing._id += 1
        self.dtype = type("DualRCGraphRingElement", (DualRCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtwtrwystoa", self._ID))

    def polynomial_coaction(self, elem: RCGraph, left=True):
        result = PA.zero @ self.zero
        for i in range(len(elem)):
            pair = elem.vertical_cut(i)
            if left:
                result += PA(*pair[0].length_vector) @ self(pair[1])
            else:
                result += self(pair[0]) @ PA(*pair[1].length_vector)
        return result

    def rc_ring_left_act(self, act_ring_elem: RCGraphRingElement, ring_elem: RCGraph):
        # all chopping off of act_elem
        ret = self.zero
        for act_elem, coeff in act_ring_elem.items():
            rc1, rc2 = ring_elem.vertical_cut(len(act_elem))
            if rc1 == act_elem:
                ret += coeff * self(rc2)
        return ret

    def rc_ring_right_act(self, ring_elem: RCGraph, act_ring_elem: RCGraphRingElement):
        # all chopping off of act_elem
        ret = self.zero
        for act_elem, coeff in act_ring_elem.items():
            rc1, rc2 = ring_elem.vertical_cut(len(act_elem))
            if rc2 == act_elem:
                ret += coeff * self(rc1)
        return ret

    def elem_sym(self, descent, weight):
        from schubmult import uncode

        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(uncode([0] * (descent - sum(weight)) + [1] * sum(weight)), len(weight), weight=weight), 1))

    def full_elem_sym(self, degree, descent, length):
        from schubmult import uncode

        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(uncode([0] * (descent - degree) + [1] * degree), length), 1))

    @property
    def zero_monom(self):
        return RCGraph([])

    def from_polynomial_algebra_element(self, elem):
        schub_elem = elem.change_basis(SchubertPolyBasis(elem.ring.genset))
        result = self.zero
        for (perm, length), coeff in schub_elem.items():
            result += coeff * self.schub(perm, length)
        return result

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

    def mul_pair(self, u_rc, v_rc):
        if u_rc.perm.inv == 0:
            return self(v_rc)
        if v_rc.perm.inv == 0:
            return self(u_rc)
        assert len(u_rc) == len(v_rc), "Multiplication only defined for RC graphs of the same number of rows"
        if u_rc.is_full_grass and v_rc.is_full_grass:
            return self(u_rc.squash_product(v_rc))
        if u_rc.is_full_grass:
            return self(v_rc.left_squash(u_rc))
        if v_rc.is_full_grass:
            return self(u_rc.squash_product(v_rc))
        u_rc_base, u_rc_grass = u_rc.squash_decomp()
        if u_rc_grass.perm.inv == 0:
            v_rc_base, v_rc_grass = v_rc.squash_decomp()
            return (self(u_rc_base.resize(len(u_rc_base) - 1)) * self(v_rc_base.resize(len(v_rc_base) - 1))).resize(len(u_rc)) * self(v_rc_grass)
        right_factor = self.mul_pair(u_rc_grass, v_rc)
        return self.mul_pair(u_rc_base, next(iter(right_factor)))
        # from .rc_graph_ring import RCGraphRing
        # if len(u_rc) != len(v_rc):
        #     return self.zero
        # if u_rc.perm.inv == 0:
        #     return self(v_rc)
        # if v_rc.perm.inv == 0:
        #     return self(u_rc)
        # if len(u_rc) == 2:
        #     u_base, u_grass = u_rc.squash_decomp()
        #     v_base, v_grass = v_rc.squash_decomp()
        #     if u_base.perm.inv == 0 and v_base.perm.inv == 0:
        #         return self(u_grass.squash_product(v_grass))
        #     if u_base.perm.inv == 0:
        #         return self(v_rc.left_squash(u_grass))
        #     if v_base.perm.inv == 0:
        #         return self(u_rc.squash_product(v_grass))
        #     # if u_grass.perm.inv == 0:
        #     #     return RCGraph([(2, 1), ()]).squash_product(v_grass)

        #     middle = v_base.left_squash(u_grass)
        #     middle_base, middle_grass = middle.squash_decomp()
        #     if middle_base.perm.inv == 0:
        #         return self(u_base.squash_product(middle_grass.squash_product(v_grass)))
        #     return self(RCGraph([(2, 1), ()]).squash_product(middle_grass.squash_product(v_grass)))
        # r = RCGraphRing()
        # return self(r(u_rc) % r(v_rc))


    def mul(self, elem1, elem2):
        if isinstance(elem2, DualRCGraphRingElement):
            result = self.zero
            for a, coeff_a in elem1.items():
                for b, coeff_b in elem2.items():

                    result += coeff_a * coeff_b * self.mul_pair(a, b)

            return result
        try:
            result_dict = {k: v * elem2 for k, v in elem1.items()}
            return self.from_dict(result_dict)
        except Exception:
            raise NotImplementedError("Multiplication with non-DualRCGraphRingElement not implemented")

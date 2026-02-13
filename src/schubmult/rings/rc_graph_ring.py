from functools import cache

from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.rings.schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.symbolic import S

from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class RCGraphRingElement(CrystalGraphRingElement, SchubertMonomialRingElement):
    """
    RCGraphRing elements are linear combinations of RCGraph basis elements.

    The product % is the polynomial product. Currently only defined when the right side
    is a dominant RC graph.

    The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
    to do this directly.

    The product * is well defined for any pair of RC graphs and is the dual product.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------

    @property
    def vex(self):
        ret = self.ring.zero
        for rc, coeff in self.items():
            ret += coeff * self.ring(rc.vex)
        return ret

    @property
    def grass(self):
        ret = self.ring.zero
        for rc, coeff in self.items():
            ret += coeff * self.ring(rc.grass)
        return ret

    def __mod__(self, other):
        """
        Polynomial product: self % other.
        Currently only defined when `other` is a dominant RC graph.
        """
        return self.ring.rc_product(self, other)

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

    # ----------------------
    # Coalgebra helpers (already present)
    # ----------------------
    def vertical_coproduct(self):
        """
        Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
        and induces the mul product.
        """
        tring = self.ring @ self.ring
        res = tring.zero
        for rc_graph, coeff in self.items():
            for i in range(len(rc_graph) + 1):
                rc1, rc2 = rc_graph.vertical_cut(i)
                res += tring.from_dict({(rc1, rc2): coeff})
        return res

    def coproduct(self):
        tring = self.ring @ self.ring
        res = tring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring.coproduct_on_basis(rc_graph)
        return res

    # ----------------------
    # Linearized crystal operations
    # ----------------------

    def raising_operator(self, index):
        """
        Linear extension of RCGraph.raising_operator:
        Apply raising_operator(index) to every basis RCGraph in self, collect results.
        Returns an RCGraphRingElement (possibly zero).
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc = rc_graph.raising_operator(index)
            if new_rc is not None:
                res += coeff * self.ring(new_rc)
        return res

    def lowering_operator(self, index):
        """
        Linear extension of RCGraph.lowering_operator.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc = rc_graph.lowering_operator(index)
            if new_rc is not None:
                res += coeff * self.ring(new_rc)
        return res

    def phi(self, index):
        """
        phi(element) := max_{basis rc in supp(element)} phi(rc)
        If element is zero, returns 0.
        """
        if len(self) == 0:
            return 0
        m = 0
        for rc_graph in self.keys():
            try:
                v = rc_graph.phi(index)
            except Exception:
                v = 0
            if v > m:
                m = v
        return m

    def epsilon(self, index):
        """
        epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)
        """
        if len(self) == 0:
            return 0
        m = 0
        for rc_graph in self.keys():
            try:
                v = rc_graph.epsilon(index)
            except Exception:
                v = 0
            if v > m:
                m = v
        return m

    def crystal_length(self):
        """
        Use maximum crystal length of basis graphs in support (0 for the zero element).
        """
        if len(self) == 0:
            return 0
        return max(getattr(rc_graph, "crystal_length", lambda: len(rc_graph))() for rc_graph in self.keys())

    def almosteq(self, other):
        return all(v == 0 for v in (self - other).values())

    def to_highest_weight(self):
        """
        Iteratively raise the element until no further raising is possible.
        Returns (highest_weight_element, raise_seq).

        Behavior notes:
        - This is the natural linear-extension of CrystalGraph.to_highest_weight.
        - The returned `highest_weight_element` is an RCGraphRingElement.
        - raise_seq is the sequence of row indices applied (in order).
        """
        rc_elem = self.ring.zero
        for rc, coeff in self.items():
            rc_elem += coeff * self.ring(rc.to_highest_weight()[0])
        return rc_elem, None

    def to_lowest_weight(self):
        """
        Iteratively raise the element until no further raising is possible.
        Returns (highest_weight_element, raise_seq).

        Behavior notes:
        - This is the natural linear-extension of CrystalGraph.to_highest_weight.
        - The returned `highest_weight_element` is an RCGraphRingElement.
        - raise_seq is the sequence of row indices applied (in order).
        """
        rc_elem = self.ring.zero
        for rc, coeff in self.items():
            rc_elem += coeff * self.ring(rc.to_lowest_weight()[0])
        return rc_elem, None

    def reverse_raise_seq(self, raise_seq):
        """
        Apply lowering_operator in reverse order to `raise_seq`.
        If the path dies (result is zero), return None (mirrors scalar behavior).
        """
        rc_elem = self
        for row in reversed(raise_seq):
            rc_elem = rc_elem.lowering_operator(row)
            if rc_elem is None or len(rc_elem) == 0:
                return None
        return rc_elem

    def crystal_reflection(self, index):
        """
        Linear extension of RCGraph.crystal_reflection:
        For each basis RCGraph, apply its crystal_reflection(index) and collect results.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.crystal_reflection(index))
        return res

    def weight_reflection(self, index):
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.weight_reflection(index))
        return res

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


class RCGraphRing(SchubertMonomialRing, CrystalGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = RCGraphRing._id
        RCGraphRing._id += 1
        self.dtype = type("RCGraphRingElement", (RCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtystoa", self._ID))

    def monomial(self, *tup):
        elem = self.one
        for a in tup:
            elem = elem * self(RCGraph.one_row(a))
        return elem

    def elem_sym(self, descent, weight):
        from schubmult import uncode

        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(uncode([0] * (descent - sum(weight)) + [1] * sum(weight)), len(weight), weight=weight), 1))

    @property
    def zero_monom(self):
        return RCGraph([])

    def from_free_algebra_element(self, elem):
        wordelem = elem.change_basis(WordBasis)
        result = self.zero
        for word, coeff in wordelem.items():
            res = self(RCGraph([]))
            for a in reversed(word):
                res = self(RCGraph.one_row(a)) * res
            result += coeff * res
        return result

    def schub(self, perm, n=None):
        """
        Return the RCGraphRing element corresponding to the Schubert polynomial
        indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).
        """
        if n is None:
            n = len(perm.trimcode)
        return self.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, n), S.One))

    # def weight_coproduct(self, elem):

    def coproduct_on_basis(self, elem):
        from schubmult import ASx

        model = ASx(elem.perm, len(elem)).coproduct()
        T = self @ self
        res = T.zero
        for (key1, key2), coeff in model.items():
            res += coeff * self.from_free_algebra_element(ASx(*key1)) @ self.from_free_algebra_element(ASx(*key2))
        res += (self._word_rc(elem.length_vector) - self(elem)) @ self(RCGraph([]).resize(len(elem))) + self(RCGraph([]).resize(len(elem))) @ (self._word_rc(elem.length_vector) - self(elem))
        return res
        # self.from_free_algebra_element(elem.perm) - elem

    # def dual_schubert_element(self, rc):
    #     from schubmult import ASx, WordBasis
    #     model = ASx(rc.perm, len(rc)).change_basis(WordBasis)
    #     for word, coeff in model.items():
    #         if coeff != 0:
    #             return self.from_free_algebra_element(WordBasis(word))

    def old_coproduct_on_basis(self, elem):
        # if not elem.is_principal:
        #     raise NotImplementedError
        tring = self @ self
        if elem.inv == 0:
            return tring.ext_multiply(self(elem), self(elem))
        # if hw:
        #     basis_elem, raise_seq = elem.to_highest_weight()
        # else:
        basis_elem = elem
        cprod = tring.zero
        right_side = True

        if right_side:
            p = basis_elem.length_vector[-1]
        else:
            p = basis_elem.length_vector[0]

        for j in range(p + 1):
            cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        if len(basis_elem) == 1:
            return cprod

        if right_side:
            lower_graph = basis_elem.vertical_cut(len(basis_elem) - 1)[0]
        else:
            lower_graph = basis_elem.rowrange(1)

        lower_module1 = self.coproduct_on_basis(lower_graph)

        if right_side:
            ret_elem = lower_module1 * cprod
        else:
            ret_elem = cprod * lower_module1

        from schubmult import CrystalGraphTensor, FreeAlgebra, Permutation, SchubertBasis  # noqa: F401

        ASx = FreeAlgebra(SchubertBasis)  # noqa: F841

        if right_side:
            up_elem2 = self(lower_graph) * self(RCGraph.one_row(p))
        else:
            up_elem2 = self(RCGraph.one_row(p)) * self(lower_graph)

        for key, coeff in up_elem2.items():
            if key != elem:
                ret_elem += (self.dual_schubert_element(key.perm, len(key)) - self(RCGraph.principal_rc(key.perm, len(key)))) @ self(RCGraph([]).resize(len(elem))) + self(
                    RCGraph([]).resize(len(elem))
                ) @ (self(key) - self(RCGraph.principal_rc(key.perm, len(key))))
                key_coprod = self.coproduct_on_basis(RCGraph.principal_rc(key.perm, len(key)))
                ret_elem -= key_coprod
        return ret_elem

    def rc_product(self, elem1, elem2):
        res = self.zero
        for u_rc, coeff_u in elem1.items():
            for v_rc, coeff_v in elem2.items():
                res += coeff_u * coeff_v * self.rc_single_product(u_rc, v_rc)
        return res

    def rc_single_product(self, u_rc, v_rc):
        # INSERTION WEIGHT TABLEAU
        # from symengine import S

        # from schubmult import CrystalGraphTensor, Plactic, RootTableau, Sx

        # from .schubert_ring import DSx
        # from .variables import GeneratingSet
        # z = GeneratingSet("z")
        if len(u_rc) != len(v_rc):
            raise NotImplementedError("Currently only defined for RC graphs of the same length.")
        if u_rc.perm.inv == 0:
            return self(v_rc)
        if v_rc.perm.inv == 0:
            return self(u_rc)
        if len(v_rc.perm.descents()) <= 1 and len(v_rc.perm.trimcode) >= len(u_rc.perm.trimcode):
            return self(u_rc.squash_product(v_rc))

        if v_rc.perm.is_dominant:
            # dominant case
            # dom_rc = v_rc
            # tensor_hw_map = {}
            # w_hw_map = {}
            # n = len(dom_rc) + 1
            # cheat_prod = Sx(dom_rc.perm) * Sx(u_rc.perm)
            # for u_rc_crystal in u_rc.full_crystal:
            #     tensor_hw_map[u_rc_crystal] = CrystalGraphTensor(dom_rc, u_rc_crystal).to_highest_weight()[0]
            #     for w in cheat_prod:
            #         dp_ret = u_rc_crystal.dualpieri(dom_rc.perm, w)
            #         if len(dp_ret) > 0:
            #             for dualpieri_seq in dp_ret:
            #                 if dualpieri_seq[-1].perm.inv == 0:
            #                     tensor_hw = tensor_hw_map[u_rc_crystal]
            #                     wt = tensor_hw.to_lowest_weight()[0].crystal_weight
            #                     wp_rcs = [rc for rc in RCGraph.all_rc_graphs(w, n - 1, weight=wt) if rc.is_lowest_weight]
            #                     wp_rc = wp_rcs[0]
            #                     if wp_rc.to_highest_weight()[0].crystal_weight == tensor_hw.crystal_weight:
            #                         w_hw_map[tensor_hw] = wp_rc.to_highest_weight()[0]
            #                     break
            # tensor = tensor_hw_map[u_rc]
            # tensor0 = CrystalGraphTensor(dom_rc, u_rc)
            # _, raise_seq = tensor0.to_highest_weight()

            # if tensor in w_hw_map:
            #     w_rc = w_hw_map[tensor]
            #     return self(w_rc.reverse_raise_seq(raise_seq))
            # collected_rcs = {}
            # for w in cheat_prod:
            #     for w_rc in RCGraph.all_hw_rcs(w, n - 1):
            #         if w_rc not in w_hw_map.values():
            #             collected_rcs[RootTableau.from_rc_graph(w_rc).weight_tableau] = w_rc

            # tab1 = RootTableau.from_rc_graph(dom_rc).weight_tableau
            # tab2 = RootTableau.from_rc_graph(tensor.factors[1]).weight_tableau

            # total_tab = Plactic().rs_insert(*tab1.row_word, *tab2.row_word)
            # try:
            #     return self(collected_rcs[total_tab].reverse_raise_seq(raise_seq))
            # except KeyError:
            #     total_tab = Plactic().rs_insert(*tab2.row_word, *tab1.row_word)
            #     return self(collected_rcs[total_tab].reverse_raise_seq(raise_seq))
            raise NotImplementedError("Multiplication only implemented for dominant right factor or Grassmannian with large enough descent.")
        raise NotImplementedError("Multiplication only implemented for dominant right factor or Grassmannian with large enough descent.")

    @cache
    def potential_products(self, left, right, length):
        len0 = max(len(left.perm.trimcode), len(right.perm.trimcode))
        left_rc_hw = left
        right_rc_hw = right
        lrc = (~(left_rc_hw)).normalize()
        rrc = (~(right_rc_hw)).normalize()
        tprods = self(lrc) * self(rrc)
        tprodst = self.zero
        for rc1, coeff1 in tprods.items():
            pain = ~rc1
            if len(pain) < max(len0, len(pain.perm.trimcode)):
                pain = pain.resize(max(len0, len(pain.perm.trimcode)) + 5)
                while len(pain) > len0 and len(pain[-1]) == 0:
                    pain = pain.zero_out_last_row()
            if len(pain) == len0:  # and pain.length_vector == tuple([a+b for a,b in zip_longest(left_rc_hw.length_vector, right_rc_hw.length_vector, fillvalue=0)]):
                tprodst += coeff1 * self(pain.resize(length))
        return set(tprodst.keys())

    def potential_prodperms(self, left, right, length):
        return set({rc.perm for rc in self.potential_products(left, right, length)})

    @property
    def one(self):
        # Define the "one" element for RCGraphRing
        identity_graph = RCGraph()
        return self.from_dict({identity_graph: 1})


class GrassRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = GrassRCGraphRing._id
        GrassRCGraphRing._id += 1
        self.dtype = type("GrassRCGraphRingElement", (RCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtasfsfystoa", self._ID))

    def new(self, x):
        return self.from_dict({x: 1})

    def __call__(self, x):
        return self.new(x)

    @property
    def zero_monom(self):
        return RCGraph([])

    def mul(self, a, b):
        return self.from_dict(super().mul(a, b))

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b))

    def from_dict(self, dct):
        dct0 = {k: v for k, v in dct.items() if len(k.perm.descents()) <= 1 and (k.perm.inv == 0 or len(k.perm.trimcode) == len(k))}
        elem = super().from_dict(dct0)
        elem.ring = self
        return elem

    def coproduct_on_basis(self, elem):
        if len(elem) == 0:
            return self(elem) @ self(elem)
        if len(elem) == 1:
            deg = elem.perm.inv
            res = (self @ self).zero
            for p in range(deg + 1):
                res += self(RCGraph.one_row(p)) @ self(RCGraph.one_row(deg - p))
            return res
        if len(elem.perm.descents()) > 1 or (len(elem) != len(elem.perm.trimcode) and elem.perm.inv > 0):
            return self.zero
        mid = len(elem) // 2
        left_cut, right_cut = elem.vertical_cut(mid)
        left_cprod = self.coproduct_on_basis(left_cut)
        right_cprod = self.coproduct_on_basis(right_cut)
        cprod = (self @ self).zero
        candidate = left_cprod * right_cprod
        for (rc1, rc2), coeff in candidate.items():
            # Filter to only Grassmannian RCs
            if next(iter((self(rc1) % self(rc2)).keys())) == elem:
                cprod += self(rc1) @ self(rc2)
        if cprod == (self @ self).zero:
            raise ValueError(f"Failed to find coproduct for {elem}")
        return cprod

    def add(self, a, b):
        return self.from_dict(super().add(a, b))

    def sub(self, a, b):
        return self.from_dict(super().sub(a, b))

    def radd(self, a, b):
        return self.from_dict(super().radd(a, b))

    def rsub(self, a, b):
        return self.from_dict(super().rsub(a, b))

    def neg(self, a):
        return self.from_dict(super().neg(a))

    # def rc_product(self, elem1, elem2):
    @property
    def zero(self):
        return self.from_dict({})

    @property
    def one(self):
        identity_graph = RCGraph()
        return self.from_dict({identity_graph: 1})

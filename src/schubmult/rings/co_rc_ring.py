from functools import cache

from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.symbolic import S, sympy_Mul

from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class CoRCGraphPrintingTerm(TypedPrintingTerm):
    # def _sympystr(self, printer=None):
    #     return printer._print(self._key)

    # def _pretty(self, printer):
    #     return printer._print(self._key)
    pass


class CoRCGraphRingElement(CrystalGraphRingElement):
    """
    CoRCGraphRing elements are linear combinations of RCGraph basis elements.
    This class implements crystal operations by linear extension: for each
    basis RCGraph in the support apply the RCGraph operation and collect
    the results into a new ring element.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------
    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k)) for k in self.keys()]

    # ----------------------
    # Coalgebra helpers (already present)
    # ----------------------
    def vertical_coproduct(self):
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
        Returns an CoRCGraphRingElement (possibly zero).
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc = rc_graph.transpose(rc_graph.rows).raising_operator(index).transpose(len(rc_graph))
            if new_rc is not None:
                res += coeff * self.ring(new_rc)
        return res

    def lowering_operator(self, index):
        """
        Linear extension of RCGraph.lowering_operator.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc = rc_graph.transpose(rc_graph.rows).lowering_operator(index).transpose(len(rc_graph))
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
                v = rc_graph.traponse(rc_graph.rows).phi(index)
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
                v = rc_graph.transpose(rc_graph.rows).epsilon(index)
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

    def to_highest_weight(self):
        """
        Iteratively raise the element until no further raising is possible.
        Returns (highest_weight_element, raise_seq).

        Behavior notes:
        - This is the natural linear-extension of CrystalGraph.to_highest_weight.
        - The returned `highest_weight_element` is an CoRCGraphRingElement.
        - raise_seq is the sequence of row indices applied (in order).
        """
        rc_elem = self
        raise_seq = []
        found = True
        while found:
            found = False
            # iterate over possible rows: 1..crystal_length()-1 (mimic scalar version)
            for row in range(1, rc_elem.crystal_length()):
                rc0 = rc_elem.raising_operator(row)
                # treat empty/zero as no raise; require that rc0 differs from rc_elem
                if rc0 is not None and len(rc0) > 0 and rc0 != rc_elem:
                    found = True
                    rc_elem = rc0
                    raise_seq.append(row)
                    break
        return rc_elem, tuple(raise_seq)

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
        For each basis CoRCGraph, apply its crystal_reflection(index) and collect results.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.crystal_reflection(index))
        return res

    def polyvalue(self, genset):
        from symengine import S
        res = S.Zero
        for rc, coeff in self.items():
            res += coeff * rc.polyvalue(genset)
        return res

# class CoRCGraphRingElement(BaseSchubertElement, CrystalGraph):
#     # def as_expr(self):
#     #     from sympy import Add

#     #     terms = []
#     #     for rc_graph, coeff in self.items():
#     #         if coeff == 1:
#     #             terms.append(self.ring.printing_term(rc_graph))
#     #         else:
#     #             terms.append(coeff * self.ring.printing_term(rc_graph))
#     #     return Add(*terms)

#     def as_ordered_terms(self, *_, **__):
#         if len(self.keys()) == 0:
#             return [S.Zero]
#         return [
#             self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k))
#             for k in self.keys()
#         ]

#     def vertical_coproduct(self):
#         tring = self.ring @ self.ring
#         res = tring.zero
#         for rc_graph, coeff in self.items():
#             for i in range(len(rc_graph) + 1):
#                 rc1, rc2 = rc_graph.vertical_cut(i)
#                 res += tring.from_dict({(rc1, rc2): coeff})
#         return res

#     def coproduct(self):
#         tring = self.ring @ self.ring
#         res = tring.zero
#         for rc_graph, coeff in self.items():
#             res += coeff * self.ring.coproduct_on_basis(rc_graph)
#         return res

#     def crystal_reflection(self, row):
#         res = self.ring.zero
#         for rc_graph, coeff in self.items():
#             res += coeff * self.ring(rc_graph.crystal_reflection(row))
#         return res


class CoRCGraphRing(CrystalGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = CoRCGraphRing._id
        CoRCGraphRing._id += 1
        self.dtype = type("CoRCGraphRingElement", (CoRCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtystoa", self._ID))

    @property
    def zero_monom(self):
        return RCGraph([])

    def printing_term(self, key):
        return CoRCGraphPrintingTerm(key)

    # def dtype(self):
    #     elem = CoRCGraphRingElement()
    #     elem.ring = self
    #     return elem

    def from_dict(self, dct):
        elem = self.dtype()
        elem.update(dct)
        return elem

    def __call__(self, key):
        return self.from_dict({key: 1})

    def from_free_algebra_element(self, elem):
        wordelem = elem.change_basis(WordBasis)
        result = self.zero
        for word, coeff in wordelem.items():
            res = self(RCGraph([]))
            for a in reversed(word):
                res = self(RCGraph.one_row(a)) * res
            result += coeff * res
        return result

    # def _one_row_cp(self, p):
    #     tring = CrystalTensorRing(self, self)
    #     res = tring.zero
    #     for i in range(p + 1):
    #         rc1, rc2 = RCGraph.one_row(i), RCGraph.one_row(p - i)
    #         res += tring.from_dict({(rc1, rc2): 1})
    #     return res

    # def _word_rc(self, word):
    #     res = self(RCGraph([]))
    #     for a in reversed(word):
    #         res *= self(RCGraph.one_row(a))
    #     return res

    # def _word_cp(self, word):
    #     tring = CrystalTensorRing(self, self)
    #     res = tring.one
    #     for a in word:
    #         res *= self._one_row_cp(a)
    #     return res

    def coproduct_on_basis(self, elem):
        if not elem.is_principal:
            raise NotImplementedError
        tring = self @ self
        if elem.inv == 0:
            return tring.ext_multiply(self(elem), self(elem))
        # if hw:
        #     basis_elem, raise_seq = elem.to_highest_weight()
        # else:
        basis_elem = elem
        cprod = tring.zero
        p = basis_elem.length_vector[-1]

        for j in range(p + 1):
            cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        if len(basis_elem) == 1:
            return cprod

        lower_graph = basis_elem.vertical_cut(len(basis_elem)-1)[0]
        lower_module1 = self.coproduct_on_basis(lower_graph)

        ret_elem = lower_module1 * cprod

        # ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if elem.perm in self.potential_prodperms(rc1, rc2, len(elem)) or elem.perm in self.potential_prodperms(rc2, rc1, len(elem))})

        # rectify
        #coprod_elem = ASx(self.perm, len(self)).coproduct()
        
        # do_it = tring.zero
        # for (rc1, rc2), _ in ret_elem.items():
        #     inner_shape = rc2.weight_tableau.shape
        #     stab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, ~rc1.perm, inner_shape)
        #     for nil_tab in stab_set:
        #         straight_tab = nil_tab.rectify()
        #         if straight_tab == rc1.p_tableau:
        #             do_it += tring((rc1, rc2))
        # ret_elem = do_it
        from schubmult import CrystalGraphTensor, Permutation
        up_elem2 =  self(lower_graph) * self(RCGraph.one_row(p))
        #elem_hw, raise_seq = elem.to_highest_weight()
        #old_ret_elem = ret_elem
        #ret_elem = tring.zero
        w0_prin = RCGraph.principal_rc(Permutation.w0(len(elem.perm)), len(elem) + 1)
        # found_one = True
        # while found_one:
        #     found_one = False
        #     for (rc1, rc2), coeff in ret_elem.items():
        #         if not CrystalGraphTensor(w0_prin, rc2.extend(1)).is_lowest_weight:
        #             ret_elem -= tring((rc1, rc2))
        #             found_one = True
        #             break
        #ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(basis_elem.perm) and k[1].perm.bruhat_leq(basis_elem.perm)})
        for key, coeff in up_elem2.items():
            if key.perm != elem.perm:
                key_coprod = self.coproduct_on_basis(RCGraph.principal_rc(key.perm, len(elem)))
                
                for (rc1_bad, rc2_bad), coeff1 in key_coprod.items():
                    #if CrystalGraphTensor(w0_prin, rc2_bad).is_lowest_weight:
                        for (rc1, rc2), coeff in ret_elem.items():
            #good_weight = tuple(Cr)
                            
                            if rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm:
            #rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm and CrystalGraphTensor(w0_prin, rc2).to_highest_weight()[0].crystal_weight != elem.length_vector:             
                                assert coeff == 1
                                if True:#CrystalGraphTensor(w0_prin, rc1_bad.extend(1)).to_lowest_weight()[0] == CrystalGraphTensor(w0_prin, rc1.extend(1)).to_lowest_weight()[0]:
                                #lw_rc2 = CrystalGraphTensor(w0_prin, rc2_bad).to_lowest_weight()[0].factors[1]
                                #lw_rc22 = CrystalGraphTensor(w0_prin, rc2).to_lowest_weight()[0].factors[1]
                                #if rc2_bad == lw_rc22:
                                    ret_elem -= tring((rc1, rc2))
                                    break
                            # if not CrystalGraphTensor(w0_prin, rc2).is_lowest_weight:
                            #     ret_elem -= tring((rc1, rc2))
                            #     break
                                #break
                        # if (rc1.edelman_greene()[0] == rc1_bad.edelman_greene()[0] and rc2.edelman_greene()[0] == rc2_bad.edelman_greene()[0]):
                        #     ret_elem -= tring((rc1, rc2))
                        #     break
                        # if rc1.p_tableau == rc1_bad.p_tableau and rc2.p_tableau == rc2_bad.p_tableau:
                        #     ret_elem -= tring((rc1, rc2))
                        #     break
                        # if rc2_bad.p_tableau == rc2.p_tableau and rc1_bad.perm == rc1.perm and rc1.p_tableau.shape == rc1_bad.p_tableau.shape:
                        #     # outer_shape = key_rc.p_tableau.shape
                        #     # inner_shape = rc2.weight_tableau.shape
                        #     # stab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, rc1.perm, inner_shape)
                        #     # for nil_tab in stab_set:
                        #     #     straight_tab = nil_tab.rectify()
                        #     #     if straight_tab == rc1.p_tableau:
                        #     #         found = True
                        #     #         break
                        #     # if found:
                        #     ret_elem -= tring((rc1, rc2))
                        #     break
                        # if rc2_bad.perm == rc2.perm and rc2_bad.p_tableau.shape == rc2.p_tableau.shape and rc1_bad.perm == rc1.perm:
                        #     # outer_shape = key_rc.p_tableau.shape
                        #     # inner_shape = rc2.weight_tableau.shape
                        #     # stab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, rc1.perm, inner_shape)
                        #     # for nil_tab in stab_set:
                        #     #     straight_tab = nil_tab.rectify()
                        #     #     if straight_tab == rc1.p_tableau:
                        #     #         found = True
                        #     #         break
                        #     # if found:
                        #     ret_elem -= tring((rc1, rc2))
                        #     break
                        # if rc2_bad.pe
                        # if rc2_bad.perm == rc2.perm and rc1_bad.perm == rc1.perm:
                        #     # outer_shape = key_rc.p_tableau.shape
                        #     # inner_shape = rc2.weight_tableau.shape
                        #     # stab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, rc1.perm, inner_shape)
                        #     # for nil_tab in stab_set:
                        #     #     straight_tab = nil_tab.rectify()
                        #     #     if straight_tab == rc1.p_tableau:
                        #     #         found = True
                        #     #         break
                        #     # if found:
                        #     ret_elem -= tring((rc1, rc2))
                        #     break
                        # if rc1_bad.p_tableau == rc1.p_tableau and rc2_bad.perm == rc2.perm:
                        #     outer_shape = key_rc.p_tableau.shape
                        #     inner_shape = rc1.weight_tableau.shape
                        #     stab_set = NilPlactic.all_skew_ed_tableaux(outer_shape, rc2.perm, inner_shape)
                        #     for nil_tab in stab_set:
                        #         straight_tab = nil_tab.rectify()
                        #         if straight_tab == rc2.p_tableau:
                        #             found = True
                        #             break
                        #     if found:
                        #         ret_elem -= tring((rc1, rc2))
                        #         break
                        
                        #     break
        ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(basis_elem.perm) and k[1].perm.bruhat_leq(basis_elem.perm)})
        #assert all(C)
        return ret_elem

    # def coproduct_on_basis(self, elem):
    #     tring = RestrictedCoRCGraphTensorRing(self, self)
    #     # trivial principal case
    #     ret_elem = tring.zero
    #     elem2 = ~elem

    #     cprd = self(elem2).vertical_coproduct()

    #     for (rc1, rc2), coeff in cprd.items():
    #         rc1_1 = ~rc1
    #         if len(rc1_1) < len(elem):
    #             rc1_1 = rc1_1.extend(len(elem) - len(rc1_1))
    #         rc2_1 = ~rc2
    #         if len(rc2_1) < len(elem):
    #             rc2_1 = rc2_1.extend(len(elem) - len(rc2_1))
    #         ret_elem += coeff * tring((rc1_1, rc2_1))
    #     return ret_elem

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

    def mul(self, a, b):
        # a, b are CoRCGraphRingElemen
        if isinstance(b, CoRCGraphRingElement):
            result_dict = {}
            for g1, c1 in a.items():
                for g2, c2 in b.items():
                    if len(g1) != len(g2):
                        continue
                    # RCGraph.prod_with_rc returns a dict {CoRCGraph: coeff}
                    prod = g1.transpose().normalize().prod_with_rc(g2.transpose().normalize())
                    for g3, c3 in prod.items():
                        # while len(g3) > max(len(g1.perm.trimcode), len(g2.perm.trimcode), 1) and len(g3[-1]) == 0:
                        #     g3 = g3.zero_out_last_row()
                        g3t = g3.transpose()
                        while len(g3t) > max(len(g1.perm.trimcode),len(g2.perm.trimcode), 1) and len(g3t[-1]) == 0:
                            g3t = g3t.zero_out_last_row()
                        if len(g3t) > len(g1):
                            continue
                        if len(g3t) < len(g1):
                            g3t = g3t.resize(len(g1))
                        result_dict[g3t] = result_dict.get(g3t, 0) + c1 * c2 * c3
            # result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        # Define the "one" element for CoRCGraphRing
        identity_graph = RCGraph()
        return self.from_dict({identity_graph: 1})


# class RestrictedCoRCGraphTensorRing(CrystalTensorRing):

#     def __init__(self, r1, r2):
#         super().__init__(r1, r2)

#     def new(self, *args):
#         try:
#             rc1, rc2 = args
#         except Exception:
#             raise ValueError("RestrictedCoRCGraphTensorRing.new expects two RCGraph arguments")
#         if len(rc1) != len(rc2):
#             raise ValueError("RestrictedCoRCGraphTensorRing.new expects CoRCGraphs of the same length")
#         return self.from_dict({(rc1, rc2): 1})

#     _z_cache = {}

#     def mul(self, elem1, elem2):sssss/
#         # elem1, elem2 are RestrictedCoRCGraphTensorRing elements
#         from schubmult.schub_lib.perm_lib import Permutation
#         def right_zero_act(rc1, rc2):
#             if (rc1, rc2) in RestrictedCoRCGraphTensorRing._z_cache:
#                 return RestrictedCoRCGraphTensorRing._z_cache[(rc1, rc2)]
#             up_perms = (ASx@ASx)(((rc1.perm, len(rc1)),(rc2.perm, len(rc2)))) * (ASx@ASx)(((Permutation([]),1),(Permutation([]),1)))

#             rc_set = set()

#             for ((perm1, _),(perm2, _)), _ in up_perms.items():
#                 for rc01 in RCGraph.all_rc_graphs(perm1, len(rc1) + 1, weight=(*rc1.length_vector, 0)):
#                     if rc01.zero_out_last_row() == rc1:
#                         for rc02 in RCGraph.all_rc_graphs(perm2, len(rc2) + 1, weight=(*rc2.length_vector, 0)):
#                             if rc02.zero_out_last_row() == rc2:
#                                 rc_set.add((rc01, rc02))

#             RestrictedCoRCGraphTensorRing._z_cache[(rc1, rc2)] = rc_set
#             return rc_set
#         result = self.zero
#         for (rc1a, rc2a), c1 in elem1.items():
#             for (rc1b, rc2b), c2 in elem2.items():
#                 num_zeros = max(len(rc1b), len(rc2b), len(rc1b.perm), len(rc2b.perm))
#                 base_pair = self.new(rc1a, rc2a)
#                 buildup = c1*c2*base_pair

#                 for _ in range(num_zeros):
#                     new_buildup = self.zero
#                     for (rc01, rc02), coeff in buildup.items():
#                         new_buildup += self.from_dict(dict.fromkeys(right_zero_act(rc01, rc02), coeff))
#                     buildup = new_buildup

#                 for (rc01, rc02), coeff in buildup.items():
#                     new_rc1 = RCGraph([*rc01[: len(rc1a)], *rc1b.shiftup(len(rc1a))])
#                     new_rc2 = RCGraph([*rc02[: len(rc2a)], *rc2b.shiftup(len(rc2a))])
#                     if new_rc1.is_valid and len(new_rc1.perm.trimcode) <= len(new_rc1) and new_rc2.is_valid and len(new_rc2.perm.trimcode) <= len(new_rc2):
#                         result += coeff * self.new(new_rc1, new_rc2)

#         return result

# def tensor_to_highest_weight(tensor_elem):
#     ret_elem = tensor_elem.ring.zero
#     for (g1, g2), coeff in tensor_elem.items():
#         (g1_hw, g2_hw), seq2 = RCGraph.to_highest_weight_pair(g1, g2)
#         ret_elem += coeff * tensor_elem.ring((g1_hw, g2_hw))
#     return ret_elem

# def tensor_to_highest_weight2(tensor_elem):
#     ret_elem = tensor_elem.ring.zero
#     for (g1, g2), coeff in tensor_elem.items():
#         g1_hw = RCGraph.principal_rc(g1.perm, len(g1))
#         g2_hw = RCGraph.principal_rc(g2.perm, len(g2))
#         ret_elem += coeff * tensor_elem.ring((g1_hw, g2_hw))
#     return ret_elem

# def ring_elem_to_highest_weight(ring_elem):
#     ret_elem = ring_elem.ring.zero
#     for g, coeff in ring_elem.items():
#         g_hw, seq = g.to_highest_weight()
#         ret_elem += coeff * ring_elem.ring(g_hw)
#     return ret_elem

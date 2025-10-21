from functools import cache

from schubmult import ASx
from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.rings.rc_graph import RCGraph
from schubmult.symbolic import S, sympy_Mul

from .crystal_graph import CrystalGraph
from .crystal_tensor_ring import CrystalTensorRing

#weight wt
# yw highest weight
# u # yv
# yv highest weight


class RCGraphPrintingTerm(TypedPrintingTerm):
    # def _sympystr(self, printer=None):
    #     return printer._print(self._key)

    # def _pretty(self, printer):
    #     return printer._print(self._key)
    pass

class RCGraphRingElement(BaseSchubertElement, CrystalGraph):
    """
    RCGraphRing elements are linear combinations of RCGraph basis elements.
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
        return [
            self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k))
            for k in self.keys()
        ]

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
        tring = CrystalTensorRing(self.ring, self.ring)
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
        return max((getattr(rc_graph, "crystal_length", lambda: len(rc_graph))() for rc_graph in self.keys()))

    def to_highest_weight(self):
        """
        Iteratively raise the element until no further raising is possible.
        Returns (highest_weight_element, raise_seq).

        Behavior notes:
        - This is the natural linear-extension of CrystalGraph.to_highest_weight.
        - The returned `highest_weight_element` is an RCGraphRingElement.
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
        For each basis RCGraph, apply its crystal_reflection(index) and collect results.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            res += coeff * self.ring(rc_graph.crystal_reflection(index))
        return res

# class RCGraphRingElement(BaseSchubertElement, CrystalGraph):
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

class RCGraphRing(BaseSchubertRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = RCGraphRing._id
        RCGraphRing._id += 1
        self.dtype = type("RCGraphRingElement", (RCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtystoa", self._ID))

    @property
    def zero_monom(self):
        return RCGraph([])

    def printing_term(self, key):
        return RCGraphPrintingTerm(key)

    # def dtype(self):
    #     elem = RCGraphRingElement()
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

    def _one_row_cp(self, p):
        tring = CrystalTensorRing(self, self)
        res = tring.zero
        for i in range(p + 1):
            rc1, rc2 = RCGraph.one_row(i), RCGraph.one_row(p - i)
            res += tring.from_dict({(rc1, rc2): 1})
        return res

    def coproduct_on_basis(self, elem):
        # we have two coproducts: weight and schub
        from schubmult.rings.crystal_graph import CrystalGraphTensor
        tring = CrystalTensorRing(self, self)

        a_coprod = ASx(elem.perm, len(elem)).coproduct()

        w_coprod = tring.one

        elem_hw, raise_seq = elem.to_highest_weight()

        for letter in elem_hw.length_vector:
            w_coprod *= self._one_row_cp(letter)
        result = tring.zero
        reserve_coprod = tring.one

        for letter in elem.length_vector:
            reserve_coprod *= self._one_row_cp(letter)
        cfs = {}
        for tensor in sorted(w_coprod.keys()):
            key = ((tensor[0].perm, len(tensor[0])), (tensor[1].perm, len(tensor[1])))
            if key not in cfs:
                cfs[key] = 0
            if key in a_coprod and cfs[key] < a_coprod[key]:
                result += tring(tensor)
                cfs[key] += 1
        lower_result = tring.zero
        for tensor, coeff in result.items():
            assert coeff == 1
            tensor = CrystalGraphTensor(*tensor)
            lowered_tensor = tensor.reverse_raise_seq(raise_seq)
            if lowered_tensor is not None:
                lower_result += coeff * tring.from_dict({lowered_tensor.factors: 1})
            else:
                last_resort = CrystalGraphTensor(*min(tensor2 for tensor2 in reserve_coprod.keys() if tensor2[0].perm == tensor.factors[0].perm and tensor2[1].perm == tensor.factors[1].perm))
                reserve_coprod -= tring.from_dict({last_resort.factors: 1})
                lower_result += coeff * tring.from_dict({last_resort.factors: 1})
        return lower_result

    def old_coproduct_on_basis(self, elem):
        tring = RestrictedRCGraphTensorRing(self, self)
        # trivial principal case
        if elem.perm.inv == 0:
            return tring((elem, elem))
        cprod = tring.zero

        p = elem.length_vector[0]
        for j in range(p + 1):
            cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        if len(elem) == 1:
            return cprod

        hw_elem, raise_seq = elem.to_lowest_weight()
        if hw_elem != elem:
            hw_coprod = self.old_coproduct_on_basis(hw_elem)
            return hw_coprod.reverse_lower_seq(raise_seq)
        # if we can multiply by zero we rule the world

        lower_elem = elem.rowrange(1, len(elem))
        lower_coprod = self.old_coproduct_on_basis(lower_elem)
        # # commuting h_i's
        # # cycles = first_row.perm.get_cycles()
        # # h = []
        # # for c in cycles:
        # #     h.append((len(c) - 1,max(c)))
        # # what is an h coproduct?
        # a_coprod_first_row = ASx(first_row.perm, len(elem)).coproduct()
        # first_row_coprod = tring.zero
        # for ((perm1, _), (perm2, _)), coeff in a_coprod_first_row.items():
        #     fr1 = RCGraph.all_rc_graphs(perm1, len(elem), weight=(perm1.inv, *first_row.length_vector[1:]))
        #     fr2 = RCGraph.all_rc_graphs(perm2, len(first_row), weight=(perm2.inv, *first_row.length_vector[1:]))
        #     first_row_coprod += coeff * tring((fr1, fr2))
        # return first_row_coprod * upper_coprod # WRONG

        candidate_coprod = cprod * lower_coprod

        # this will have extra terms we need to remove
        ret_elem = candidate_coprod

        if p == 0:
            return ret_elem

        up_a_elem = self(RCGraph.one_row(p)) * self(lower_elem)

        if len(up_a_elem) == 1:
            return ret_elem

        for key_rc, coeff in up_a_elem.items():
            if key_rc.perm != elem.perm:
                assert coeff == 1
                # a_elem = ASx(key_rc.perm, len(key_rc)).coproduct()
                # key_rc_hw, raise_seq = key_rc.to_highest_weight()
                # if key_rc_hw != key_rc:
                #     print("pink")
                #     bad_coprod = self(key_rc_hw).coproduct().reverse_raise_seq(raise_seq)
                #     ret_elem -= bad_coprod
                # else:
                bad_coprod = self.old_coproduct_on_basis(RCGraph.principal_rc(key_rc.perm,len(key_rc)))
                #bad_coprod = self.coproduct_on_basis(key_rc)
                # ret_elem -= bad_coprod
                for (rc1_bad, rc2_bad), cff2 in bad_coprod.items():
                    ret_elem -= cff2 * tring.new(*min([(rc1, rc2) for (rc1, rc2) in ret_elem.keys() if rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm]))

        # test_elem = tring.zero
        # for (rc1, rc2), v in ret_elem.items():
        #     test_elem += v * tring((_trim_and_zero(rc1), _trim_and_zero(rc2)))
        # assert test_elem == lower_elem, f"Coproduct check failed for {elem=}: got {test_elem=}, expected {lower_elem=}"

        ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(elem.perm) and rc2.perm.bruhat_leq(elem.perm)})

        return ret_elem

    def mul(self, a, b):
        # a, b are RCGraphRingElemen
        if isinstance(b, RCGraphRingElement):
            result_dict = {}
            for g1, c1 in a.items():
                for g2, c2 in b.items():
                    # RCGraph.prod_with_rc returns a dict {RCGraph: coeff}
                    prod = g1.prod_with_rc(g2)
                    for g3, c3 in prod.items():
                        result_dict[g3] = result_dict.get(g3, 0) + c1 * c2 * c3
            # result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        # Define the "one" element for RCGraphRing
        identity_graph = RCGraph()
        return self.from_dict({identity_graph: 1})

class RestrictedRCGraphTensorRing(CrystalTensorRing):

    def __init__(self, r1, r2):
        super().__init__(r1, r2)

    def new(self, *args):
        try:
            rc1, rc2 = args
        except Exception:
            raise ValueError("RestrictedRCGraphTensorRing.new expects two RCGraph arguments")
        if len(rc1) != len(rc2):
            raise ValueError("RestrictedRCGraphTensorRing.new expects RCGraphs of the same length")
        return self.from_dict({(rc1, rc2): 1})

    _z_cache = {}

    def mul(self, elem1, elem2):
        # elem1, elem2 are RestrictedRCGraphTensorRing elements
        from schubmult.perm_lib import Permutation
        def right_zero_act(rc1, rc2):
            if (rc1, rc2) in RestrictedRCGraphTensorRing._z_cache:
                return RestrictedRCGraphTensorRing._z_cache[(rc1, rc2)]
            up_perms = (ASx@ASx)(((rc1.perm, len(rc1)),(rc2.perm, len(rc2)))) * (ASx@ASx)(((Permutation([]),1),(Permutation([]),1)))

            rc_set = set()

            for ((perm1, _),(perm2, _)), _ in up_perms.items():
                for rc01 in RCGraph.all_rc_graphs(perm1, len(rc1) + 1, weight=(*rc1.length_vector, 0)):
                    if rc01.zero_out_last_row() == rc1:
                        for rc02 in RCGraph.all_rc_graphs(perm2, len(rc2) + 1, weight=(*rc2.length_vector, 0)):
                            if rc02.zero_out_last_row() == rc2:
                                rc_set.add((rc01, rc02))

            RestrictedRCGraphTensorRing._z_cache[(rc1, rc2)] = rc_set
            return rc_set
        result = self.zero
        for (rc1a, rc2a), c1 in elem1.items():
            for (rc1b, rc2b), c2 in elem2.items():
                num_zeros = max(len(rc1b), len(rc2b), len(rc1b.perm), len(rc2b.perm))
                base_pair = self.new(rc1a, rc2a)
                buildup = c1*c2*base_pair

                for _ in range(num_zeros):
                    new_buildup = self.zero
                    for (rc01, rc02), coeff in buildup.items():
                        new_buildup += self.from_dict(dict.fromkeys(right_zero_act(rc01, rc02), coeff))
                    buildup = new_buildup

                for (rc01, rc02), coeff in buildup.items():
                    new_rc1 = RCGraph([*rc01[: len(rc1a)], *rc1b.shiftup(len(rc1a))])
                    new_rc2 = RCGraph([*rc02[: len(rc2a)], *rc2b.shiftup(len(rc2a))])
                    if new_rc1.is_valid and len(new_rc1.perm.trimcode) <= len(new_rc1) and new_rc2.is_valid and len(new_rc2.perm.trimcode) <= len(new_rc2):
                        result += coeff * self.new(new_rc1, new_rc2)

        return result

def tensor_to_highest_weight(tensor_elem):
    ret_elem = tensor_elem.ring.zero
    for (g1, g2), coeff in tensor_elem.items():
        (g1_hw, g2_hw), seq2 = RCGraph.to_highest_weight_pair(g1, g2)
        ret_elem += coeff * tensor_elem.ring((g1_hw, g2_hw))
    return ret_elem

def tensor_to_highest_weight2(tensor_elem):
    ret_elem = tensor_elem.ring.zero
    for (g1, g2), coeff in tensor_elem.items():
        g1_hw = RCGraph.principal_rc(g1.perm, len(g1))
        g2_hw = RCGraph.principal_rc(g2.perm, len(g2))
        ret_elem += coeff * tensor_elem.ring((g1_hw, g2_hw))
    return ret_elem

def ring_elem_to_highest_weight(ring_elem):
    ret_elem = ring_elem.ring.zero
    for g, coeff in ring_elem.items():
        g_hw, seq = g.to_highest_weight()
        ret_elem += coeff * ring_elem.ring(g_hw)
    return ret_elem



########MOST RECENT###########

# raise to highest weight and remember the raise sequence
        # basis_elem, raise_seq = elem.to_highest_weight()

        # # --- LR / crystal matching using ASx(...).coproduct() ---
        # try:
        #     fa_coprod = ASx(basis_elem.perm, len(basis_elem)).coproduct()
        # except Exception:
        #     fa_coprod = {}

        # # normalize to mapping (perm_u, perm_v) -> multiplicity
        # schubert_pairs = {}
        # for k, mult in fa_coprod.items():
        #     try:
        #         (perm_u, _), (perm_v, _) = k
        #     except Exception:
        #         continue
        #     schubert_pairs[(perm_u, perm_v)] = schubert_pairs.get((perm_u, perm_v), 0) + mult

        # target_lv = tuple(basis_elem.length_vector)
        # target_len = len(target_lv)

        # def _is_highest_weight(rcg):
        #     for r in range(1, rcg.crystal_length()):
        #         try:
        #             if rcg.raising_operator(r) is not None:
        #                 return False
        #         except Exception:
        #             return False
        #     return True

        # phi_indices = range(1, target_len)
        # #hw_tensor = tring.zero

        # for (perm_u, perm_v), mult in schubert_pairs.items():
        #     # right candidates must be highest-weight and have exact length-vector length
        #     try:
        #         right_candidates = [
        #             R for R in RCGraph.all_rc_graphs(perm_v, len(basis_elem))
        #             if _is_highest_weight(R)
        #         ]
        #     except Exception:
        #         right_candidates = []

        #     if not right_candidates:
        #         continue

        #     # left candidates with exact length-vector length
        #     try:
        #         left_candidates = [
        #             L for L in RCGraph.all_rc_graphs(perm_u, len(basis_elem))
        #             if len(tuple(L.length_vector)) == target_len
        #         ]
        #     except Exception:
        #         left_candidates = []

        #     # cheap prefilter: left length_vector entries must not exceed target entries
        #     left_candidates = [
        #         L for L in left_candidates
        #         if all(lv <= tv for lv, tv in zip(tuple(L.length_vector), target_lv))
        #     ]

        #     matching_pairs = []
        #     for R in right_candidates:
        #         # precompute phi(R) for indices of interest
        #         phi_R = {}
        #         for i in phi_indices:
        #             try:
        #                 phi_R[i] = R.phi(i)
        #             except Exception:
        #                 phi_R[i] = 0

        #         for L in left_candidates:
        #             lvL = tuple(L.length_vector)
        #             lvR = tuple(R.length_vector)

        #             # strict same-length requirement; require elementwise sum equals target
        #             if tuple(a + b for a, b in zip(lvL, lvR)) != target_lv:
        #                 continue

        #             # epsilon(L)_i <= phi(R)_i for all i in phi_indices
        #             ok = True
        #             for i in phi_indices:
        #                 eps_L = L.epsilon(i)
        #                 if eps_L > phi_R.get(i, 0):
        #                     ok = False
        #                     break
        #             if ok:
        #                 matching_pairs.append((L, R))
        # assert len(matching_pairs)>0, "matching_pairs should be defined"

        # #basis_elem, raise_seq = elem.to_highest_weight()

        # cprod = tring.zero
        # p = basis_elem.length_vector[0]

        # for j in range(p + 1):
        #     cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        # if len(basis_elem) == 1:
        #     return cprod

        # lower_graph = basis_elem.rowrange(1, len(basis_elem))
        # lower_module1 = self.coproduct_on_basis(lower_graph)

        # ret_elem = cprod * lower_module1

        # up_elem2 = self(RCGraph.one_row(p)) * self(lower_graph)
        # for key, coeff in up_elem2.items():
        #     if key.perm != basis_elem.perm:
        #         assert coeff == 1
        #         key_rc = RCGraph.principal_rc(key.perm, len(key))
        #         cp = self.coproduct_on_basis(key_rc)
        #         for (rc1_bad, rc2_bad), cff2 in cp.items():
        #             for (rc1, rc2), v in ret_elem.items():
        #                 if (rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm) and (rc1, rc2) not in matching_pairs:
        #                     ret_elem -= tring((rc1, rc2))
        #                     break


        # ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})
        # ret_elem = ret_elem.reverse_raise_seq(raise_seq)
        # # for (rc1, rc2), v in ret_elem.items():
        # #     if (rc1, rc2) not in matching_pairs:
        # #         print(f"Found non-matching pair in coproduct of {elem=}: {(rc1, rc2)=}")
        # # print(f"{matching_pairs=}")

        # return ret_elem


#################################

# @cache
    # def coproduct_on_basis(self, elem):
        # ## RIGHT ##
        # # simulate principal
        # # left = False
        # tring = self@self
        # #basis_elem, raise_seq = elem.to_highest_weight()
        # if elem.perm.inv == 0:
        #     return tring((elem, elem))
        # # print(f"Trying {elem=}")
        # basis_elem, raise_seq = elem.to_highest_weight()
        # #basis_elem = elem 
        # # print(f"Highest weight {basis_elem=}")
        # # if not basis_elem.is_principal:
        # #     raise NotImplementedError(f"Coproduct only implemented for principal highest weight rc graphs, got {basis_elem=}")

        # #basis_elem, raise_seq = basis_elem.to_highest_weight()
        # # if basis_elem != basis_elem:
        # #     print(f"We have {basis_elem=} {raise_seq=} {basis_elem=}")
        # #     assert basis_elem.length_vector != basis_elem.length_vector
        # #     # first do the coproduct on the highest weight vector
        # #     high_cop = self.coproduct_on_basis(basis_elem)
        # #     # then raise each component


        # #if basis_elem == elem:
        # # if left:
        # cprod = tring.zero
        # p = basis_elem.length_vector[0]

        # for j in range(p + 1):
        #     cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        # if len(basis_elem) == 1:
        #     return cprod

        # lower_graph = basis_elem.rowrange(1,len(basis_elem))
        # lower_module1 = self.coproduct_on_basis(lower_graph)

        # ret_elem = cprod * lower_module1

        # # ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})
        # up_elem2 = self(RCGraph.one_row(p)) * self(lower_graph)
        # #     for key, coeff in up_elem2.items():
        # #         if key.perm != basis_elem.perm:
        # #             assert coeff == 1
        # #             key_rc = RCGraph.principal_rc(key.perm, len(key))
        # #             cp = self.coproduct_on_basis(key_rc)
        # #             for (rc1_bad, rc2_bad), cff2 in cp.items():
        # #                 for (rc1, rc2), v in ret_elem.items():
        # #                     if rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm:
        # #                         ret_elem -= tring((rc1, rc2))
        # #                         break
        # # else:
        # # cprod = tring.zero
        # # p = basis_elem.length_vector[-1]

        # # for j in range(p + 1):
        # #     cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        # # if len(basis_elem) == 1:
        # #     return cprod

        # # lower_graph = basis_elem.vertical_cut(len(basis_elem) - 1)[0]
        # # lower_module1 = self.coproduct_on_basis(lower_graph)

        # # ret_elem = lower_module1 * cprod

        # ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})

        # # @cache
        # # def co_principal_perms(perm, length):
        # #     from schubmult import ASx, uncode

        # #     if len(perm.trimcode) == 0:
        # #         return []
        # #     if length < len(perm.trimcode):
        # #         return []
        # #     lower_elem = ASx(uncode(perm.trimcode[1:]), length - 1)
        # #     upper_elem = ASx(uncode([perm.trimcode[0]]), 1) * lower_elem
        # #     return [key[0] for key in upper_elem.keys() if key[0] != perm]
        # # if len(co_principal_perms(basis_elem.perm, len(basis_elem))) == 0:
        # #     return ret_elem
        # # up_elem2 = self(lower_graph) * self(RCGraph.one_row(p))
        # for key, coeff in up_elem2.items():
        #     if key.perm != basis_elem.perm:
        #         assert coeff == 1
        #         key_rc = RCGraph.principal_rc(key.perm, len(key))
        #         cp = self.coproduct_on_basis(key_rc)
        #         # ret_elem -= coeff * cp
        #         for (rc1_bad, rc2_bad), cff2 in cp.items():
        #             for (rc1, rc2), v in ret_elem.items():
        #                 if rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm:
        #                     ret_elem -= tring((rc1, rc2))
        #                     break

        # ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})

        # # else:
        # #     basis_cp = self.coproduct_on_basis(basis_elem)
        # #     ret_elem = tring.zero
        # # if len(basis_elem.perm.trimcode) == len(basis_elem):
        # #     for (rc1, rc2), cff in ret_elem.items():
        # #         (rc11, rc22), raise_seq2 = RCGraph.to_highest_weight_pair(rc1, rc2)
        # #         assert rc1 == rc11 and rc2 == rc22, f"{rc1=} {rc11=} {rc2=} {rc22=}"
        # # assert at_least_one, f"No highest weight terms found in coproduct of {elem=}, got {ret_elem=}"
        # # if basis_elem != elem:
        # #     ret_elem = tring.from_dict({(RCGraph.reverse_raise_seq_pair(rc1, rc2, raise_seq)): v for (rc1, rc2), v in ret_elem.items()})
        # #     elem2 = basis_elem.reverse_raise_seq(raise_seq)
        # #     assert elem2 == elem, f"{elem2=} {elem=}"
        # return ret_elem
        #@cache
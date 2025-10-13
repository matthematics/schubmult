from functools import cache

from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.rings.rc_graph import RCGraph


class RCGraphPrintingTerm(TypedPrintingTerm):
    def _sympystr(self, printer=None):
        return printer._print(self._key)

    def _pretty(self, printer):
        return printer._print(self._key)


class RCGraphRingElement(BaseSchubertElement):
    def as_expr(self):
        from sympy import Add

        terms = []
        for rc_graph, coeff in self.items():
            if coeff == 1:
                terms.append(self.ring.printing_term(rc_graph))
            else:
                terms.append(coeff * self.ring.printing_term(rc_graph))
        return Add(*terms)

    def vertical_coproduct(self):
        tring = self.ring @ self.ring
        res = tring.zero
        for rc_graph, coeff in self.items():
            for i in range(len(rc_graph) + 1):
                rc1, rc2 = rc_graph.vertical_cut(i)
                res += tring.from_dict({(rc1, rc2): coeff})
        return res

class RCGraphRing(BaseSchubertRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = RCGraphRing._id
        RCGraphRing._id += 1

    def __hash__(self):
        return hash(("Dinkberrtystoa", self._ID))

    @property
    def zero_monom(self):
        return RCGraph()

    def printing_term(self, key):
        return RCGraphPrintingTerm(key)

    def dtype(self):
        elem = RCGraphRingElement()
        elem.ring = self
        return elem

    def from_dict(self, dct):
        elem = self.dtype()
        elem.update(dct)
        return elem

    def __call__(self, key):
        return self.from_dict({key: 1})

    @cache
    def coproduct_on_basis(self, elem):
        ## RIGHT ##
        # simulate principal
        # left = False
        tring = self@self
        #basis_elem, raise_seq = elem.to_highest_weight()
        if elem.perm.inv == 0:
            return tring((elem, elem))
        # print(f"Trying {elem=}")
        basis_elem, raise_seq = elem.to_highest_weight()
        #basis_elem = elem 
        # print(f"Highest weight {basis_elem=}")
        # if not basis_elem.is_principal:
        #     raise NotImplementedError(f"Coproduct only implemented for principal highest weight rc graphs, got {basis_elem=}")

        #basis_elem, raise_seq = basis_elem.to_highest_weight()
        # if basis_elem != basis_elem:
        #     print(f"We have {basis_elem=} {raise_seq=} {basis_elem=}")
        #     assert basis_elem.length_vector != basis_elem.length_vector
        #     # first do the coproduct on the highest weight vector
        #     high_cop = self.coproduct_on_basis(basis_elem)
        #     # then raise each component


        #if basis_elem == elem:
        # if left:
        #     cprod = tring.zero
        #     p = basis_elem.length_vector[0]

        #     for j in range(p + 1):
        #         cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        #     if len(basis_elem) == 1:
        #         return cprod

        #     lower_graph = basis_elem.rowrange(1,len(basis_elem))
        #     lower_module1 = self.coproduct_on_basis(lower_graph)

        #     ret_elem = cprod * lower_module1

        #     ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})
        #     up_elem2 = self(RCGraph.one_row(p)) * self(lower_graph)
        #     for key, coeff in up_elem2.items():
        #         if key.perm != basis_elem.perm:
        #             assert coeff == 1
        #             key_rc = RCGraph.principal_rc(key.perm, len(key))
        #             cp = self.coproduct_on_basis(key_rc)
        #             for (rc1_bad, rc2_bad), cff2 in cp.items():
        #                 for (rc1, rc2), v in ret_elem.items():
        #                     if rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm:
        #                         ret_elem -= tring((rc1, rc2))
        #                         break
        # else:
        cprod = tring.zero
        p = basis_elem.length_vector[-1]

        for j in range(p + 1):
            cprod += tring.ext_multiply(self(RCGraph.one_row(j)), self(RCGraph.one_row(p - j)))
        if len(basis_elem) == 1:
            return cprod

        lower_graph = basis_elem.vertical_cut(len(basis_elem) - 1)[0]
        lower_module1 = self.coproduct_on_basis(lower_graph)

        ret_elem = lower_module1 * cprod

        ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})
        up_elem2 = self(lower_graph) * self(RCGraph.one_row(p))
        for key, coeff in up_elem2.items():
            if key.perm != basis_elem.perm:
                assert coeff == 1
                key_rc = RCGraph.principal_rc(key.perm, len(key))
                cp = self.coproduct_on_basis(key_rc)
                for (rc1_bad, rc2_bad), cff2 in cp.items():
                    for (rc1, rc2), v in ret_elem.items():
                        if rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm:
                            ret_elem -= tring((rc1, rc2))
                            break

        ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})

        # else:
        #     basis_cp = self.coproduct_on_basis(basis_elem)
        #     ret_elem = tring.zero
        #     for (rc1, rc2), cff in basis_cp.items():
        #         ret_elem += tring(RCGraph.reverse_raise_seq_pair(rc1, rc2, raise_seq))
        if basis_elem != elem:
            ret_elem = tring.from_dict({(RCGraph.reverse_raise_seq_pair(rc1, rc2, raise_seq)): v for (rc1, rc2), v in ret_elem.items()})
            elem2 = basis_elem.reverse_raise_seq(raise_seq)
            assert elem2 == elem, f"{elem2=} {elem=}"
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

    def zero(self):
        return self.dtype()

    def one(self):
        # Define the "one" element for RCGraphRing
        identity_graph = RCGraph()
        return self.from_dict({identity_graph: 1})

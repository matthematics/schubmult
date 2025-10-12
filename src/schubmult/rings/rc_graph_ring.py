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

    def coproduct_on_basis(self, basis_elem):
        from schubmult import ASx, uncode
        # simulate principal
        tring = self@self
        if len(basis_elem) == 0:
            return tring((RCGraph(), RCGraph()))
        if basis_elem.perm.inv == 0:
            return tring((basis_elem, basis_elem))
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
        # don't do ASx here to try
        #up_elem2 = ASx(lower_graph.perm, len(lower_graph)) * ASx(uncode([p]), 1)
        up_elem2 = self(lower_graph) * self(RCGraph.one_row(p))
        for key, coeff in up_elem2.items():
            if key.perm != basis_elem.perm:
                assert coeff == 1
                for (rc1_bad, rc2_bad), cff2 in self.coproduct_on_basis(RCGraph.principal_rc(key.perm, len(key))).items():
                #for (rc1_bad, rc2_bad), cff2 in (self.coproduct_on_basis(key.vertical_cut(len(key)-1)[0])*self.coproduct_on_basis(RCGraph.one_row(len(key[-1])))).items():
                    keys2 = set(ret_elem.keys())
                    for rc1, rc2 in keys2:
                        if (rc1.perm == rc1_bad.perm and rc2.perm == rc2_bad.perm):
                            ret_elem -= tring((rc1, rc2))
                            break
        ret_elem = tring.from_dict({(rc1, rc2): v for (rc1, rc2), v in ret_elem.items() if rc1.perm.bruhat_leq(basis_elem.perm) and rc2.perm.bruhat_leq(basis_elem.perm)})
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

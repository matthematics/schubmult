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

from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.rings.rc_graph import RCGraph


class RCGraphPrintingTerm(TypedPrintingTerm):
    pass


class RCGraphRingElement(BaseSchubertElement):
    # No need to override __mul__; handled by BaseSchubertElement
    def as_expr(self):
        from sympy import Add

        terms = []
        for rc_graph, coeff in self.items():
            if coeff == 1:
                terms.append(self.ring.printing_term(rc_graph))
            else:
                terms.append(coeff * self.ring.printing_term(rc_graph))
        return Add(*terms)


class RCGraphRing(BaseSchubertRing):
    def __init__(self, *_, **__):
        pass

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
        else:
            result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    def zero(self):
        return self.dtype()

    def one(self):
        # Define the "one" element for RCGraphRing
        identity_graph = RCGraph.identity()
        return self.dtype(self, {identity_graph: 1})

from functools import cache

from schubmult import ASx
from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.rings.plactic import NilPlactic
from schubmult.symbolic import S, sympy_Mul

#weight wt
# yw highest weight
# u # yv
# yv highest weight


class CoxeterKnuthPrintingTerm(TypedPrintingTerm):
    pass

class CoxeterKnuthRingElement(BaseSchubertElement):
    """
    CoxeterKnuthRing elements are linear combinations of NilPlactic basis elements.
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


class CoxeterKnuthRing(BaseSchubertRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = CoxeterKnuthRing._id
        CoxeterKnuthRing._id += 1
        self.dtype = type("CoxeterKnuthRingElement", (CoxeterKnuthRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtystoa", "poing", self._ID))

    @property
    def zero_monom(self):
        return (NilPlactic(), 0)

    def printing_term(self, key):
        return CoxeterKnuthPrintingTerm(key)

    # def dtype(self):
    #     elem = CoxeterKnuthRingElement()
    #     elem.ring = self
    #     return elem

    def from_dict(self, dct):
        elem = self.dtype()
        elem.update(dct)
        return elem

    def __call__(self, key):
        return self.from_dict({key: 1})

    def mul(self, a, b):
        # a, b are CoxeterKnuthRingElemen
        if isinstance(b, CoxeterKnuthRingElement):
            result_dict = {}
            for (g1, len1), c1 in a.items():
                for (g2, len2), c2 in b.items():
                    # CoxeterKnuth.prod_with_rc returns a dict {CoxeterKnuth: coeff}
                    prod = g1.hw_rc(len1).prod_with_rc(g2.hw_rc(len2))
                    for g3, c3 in prod.items():
                        result_dict[(g3.p_tableau,len(g3))] = result_dict.get((g3.p_tableau,len(g3)), 0) + c1 * c2 * c3
            # result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        # Define the "one" element for CoxeterKnuthRing
        identity_graph = (NilPlactic(), 0)
        return self.from_dict({identity_graph: 1})


from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.schub_lib.crystal_graph import CrystalGraphTensor
from schubmult.schub_lib.nilplactic import NilPlactic
from schubmult.schub_lib.plactic import Plactic
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.symbolic import S, sympy_Mul

from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement

#weight wt
# yw highest weight
# u # yv
# yv highest weight


class CoxeterKnuthPrintingTerm(TypedPrintingTerm):
    pass

class CoxeterKnuthRingElement(CrystalGraphRingElement):
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

    def __eq__(self, other):
        return type(self) is type(other) and dict(self) == dict(other)

class CoxeterKnuthKey(Plactic):

    def __init__(self, p_tableau, weight_tableau, length):
        self._p_tableau = p_tableau
        self._word = weight_tableau._word
        self._weight_tableau = weight_tableau
        self.length = length

    def crystal_length(self) -> int:
        return self.length

    def _sympystr(self, printer):
        return printer._print(CrystalGraphTensor(self._p_tableau, self._weight_tableau))

    def _pretty(self, printer):
        return printer._print(CrystalGraphTensor(self._p_tableau, self._weight_tableau))

    @property
    def rc_graph(self):
        try:
            rc = self._p_tableau.hw_rc(self.length)
        except Exception:
            return None
        if rc:
            _, raise_seq = self._weight_tableau.to_highest_weight(length=self.length)
            rc = rc.reverse_raise_seq(raise_seq)
        return rc

    @classmethod
    def from_rc_graph(cls, rc: RCGraph):
        return cls(rc.p_tableau, rc.weight_tableau, len(rc))


class CoxeterKnuthRing(CrystalGraphRing):
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
        return self.from_dict({CoxeterKnuthKey(*key): 1})

    def mul(self, a, b):
        # a, b are CoxeterKnuthRingElemen
        if isinstance(b, CoxeterKnuthRingElement):
            result_dict = {}
            for g1, c1 in a.items():
                for g2, c2 in b.items():
                    # CoxeterKnuth.product returns a dict {CoxeterKnuth: coeff}
                    prod = g1.rc_graph.product(g2.rc_graph)
                    for g3, c3 in prod.items():
                        result_dict[CoxeterKnuthKey(g3.p_tableau,g3.weight_tableau,len(g3))] = result_dict.get(CoxeterKnuthKey(g3.p_tableau,g3.weight_tableau,len(g3)), 0) + c1 * c2 * c3
            # result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    def __eq__(self, other):
        return type(self) is type(other) and self._ID == other._ID

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        # Define the "one" element for CoxeterKnuthRing
        identity_graph = (NilPlactic(), Plactic(), 0)
        return self.from_dict({CoxeterKnuthKey(*identity_graph): 1})


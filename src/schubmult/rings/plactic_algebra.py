from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from schubmult.schub_lib.plactic import Plactic
from schubmult.symbolic import S, sympy_Mul

#weight wt
# yw highest weight
# u # yv
# yv highest weight


class PlacticPrintingTerm(TypedPrintingTerm):
    pass

class PlacticAlgebraElement(CrystalGraphRingElement):
    """
    PlacticAlgebra elements are linear combinations of Plactic basis elements.
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


class PlacticAlgebra(CrystalGraphRing):
    _id = 0

    def __init__(self, *_, op=False, **__):
        self._ID = PlacticAlgebra._id
        PlacticAlgebra._id += 1
        self.dtype = type("PlacticAlgebraElement", (PlacticAlgebraElement,), {"ring": self})
        self._op = op

    def __hash__(self):
        return hash(("Dinkberrtystoa", "poing", self._ID))

    @property
    def zero_monom(self):
        return (Plactic(), 0)

    def printing_term(self, key):
        return PlacticPrintingTerm(key)

    # def dtype(self):
    #     elem = PlacticAlgebraElement()
    #     elem.ring = self
    #     return elem

    def from_dict(self, dct):
        elem = self.dtype()
        elem.update(dct)
        return elem

    def __call__(self, key):
        return self.from_dict({key: 1})

    def mul(self, a, b):
        # a, b are PlacticAlgebraElemen
        if isinstance(b, PlacticAlgebraElement):
            result_dict = {}
            if self._op:
                a, b = b, a
            for g1, c1 in a.items():
                for g2,  c2 in b.items():
                    # Plactic.prod_with_rc returns a dict {Plactic: coeff}
                    g3 = g1 * g2
                    result_dict[g3] = result_dict.get(g3, 0) + c1 * c2
            # result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    def __eq__(self, other):
        return type(self) is type(other) and self._ID == other._ID

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        # Define the "one" element for PlacticAlgebra
        identity_graph = (Plactic(), 0)
        return self.from_dict({identity_graph: 1})


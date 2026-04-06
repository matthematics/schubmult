
from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


class HWRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = HWRCGraphRing._id
        HWRCGraphRing._id += 1
        self.dtype = type("HWRCGraphRingElement", (RCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbeasfrrtasfsfystoa", self._ID))

    def new(self, x):
        return self.from_dict({x: 1})

    def __call__(self, x):
        return self.new(x)

    @property
    def zero_monom(self):
        return RCGraph([])

    def _snap_highest_weight(self, elem):
        ret = self.zero
        for key, coeff in elem.items():
            ret += self.from_dict({key.to_highest_weight()[0]: coeff})
        return ret

    def mul(self, a, b):
        return self.from_dict(self._snap_highest_weight(super().mul(a, b)))

    def rmul(self, a, b):
        return self.from_dict(self._snap_highest_weight(super().rmul(a, b)))

    def from_dict(self, dct):
        elem = super().from_dict(dct)
        elem.ring = self
        return elem

    def coproduct_on_basis(self, elem):
        from schubmult import ASx, GrassTensorAlgebra
        from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
        g = GrassTensorAlgebra()
        r = RCGraphRing()
        if elem.perm.inv == 0:
            return self(elem) @ self(elem)
        if not elem.is_highest_weight:
            raise ValueError("Can only coproduct highest weight elements")
        cprd = ASx(elem.perm, len(elem)).coproduct()
        result = (r@r).zero
        for ((perm1, _), (perm2, _)), coeff in cprd.items():
            cem1 = RCGraph.full_CEM(perm1, len(elem))
            cem2 = RCGraph.full_CEM(perm2, len(elem))
            for rc1, cem_dict1 in cem1.items():
                for rc2, cem_dict2 in cem2.items():
                    for key1, coeff1 in cem_dict1.items():
                        for key2, coeff2 in cem_dict2.items():
                            tensor = CrystalGraphTensor(*key1, *key2)
                            if tensor.is_highest_weight and self.from_dict((g(key1) * g(key2)).to_rc_graph_ring_element()).almosteq(self(elem)):
                                #result += coeff1 * coeff2 * self.from_dict(self._snap_highest_weight(g(key1).to_rc_graph_ring_element().resize(len(elem)))) @ self.from_dict(self._snap_highest_weight(g(key2).to_rc_graph_ring_element().resize(len(elem))))
                                result += coeff1 * coeff2 * g(key1).to_rc_graph_ring_element().resize(len(elem)) @ g(key2).to_rc_graph_ring_element().resize(len(elem))
        return result



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

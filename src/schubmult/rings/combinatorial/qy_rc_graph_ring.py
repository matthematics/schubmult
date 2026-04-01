
from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


class QYRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = QYRCGraphRing._id
        QYRCGraphRing._id += 1
        self.dtype = type("QYRCGraphRingElement", (RCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbeasfrrtasfasfasfsfsfsfssfystoa", self._ID))

    def new(self, x):
        return self.from_dict({x: 1})

    def __call__(self, x):
        return self.new(x)

    @property
    def zero_monom(self):
        return RCGraph([])

    def _snap_qy(self, elem):
        ret = self.zero

        def the_qy(rc):
            return next(iter([rc0 for rc0 in RCGraph.all_rc_graphs(rc.perm, len(rc)) if rc0.perm_word == rc.perm_word and rc0.is_quasi_yamanouchi]))

        for key, coeff in elem.items():
            ret += self.from_dict({the_qy(key): coeff})
        return ret

    def mul(self, a, b):
        return self.from_dict(self._snap_qy(super().mul(a, b)))

    def rmul(self, a, b):
        return self.from_dict(self._snap_qy(super().rmul(a, b)))

    def from_dict(self, dct):
        elem = super().from_dict(dct)
        elem.ring = self
        return elem


"""`QYRCGraphRing`: `RCGraphRing` quotient onto quasi-Yamanouchi RC graphs.

Every product is snapped by merging mergeable adjacent rows (``_canonical_rc``, the
same normalization as `RCGraph.snap_qy`), so basis elements are quasi-Yamanouchi.
"""

from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


def _canonical_rc(rc):
    """Merge mergeable adjacent rows until ``rc`` is quasi-Yamanouchi."""
    if rc.is_quasi_yamanouchi:
        return rc
    #row = 1
    for i in range(1, len(rc)):
        if max(rc[i], default=0) < min(rc[i - 1], default=0) and min(rc[i - 1], default=0) >= i + 1:
            new_rc = [*rc]
            new_rc[i] = (*rc[i - 1], *rc[i])
            new_rc[i - 1] = ()
            return _canonical_rc(RCGraph(new_rc))
    raise ValueError(f"Should never reach here {rc=} {rc.is_quasi_yamanouchi=}")


class QYRCGraphRing(RCGraphRing):
    """`RCGraphRing` with products snapped to quasi-Yamanouchi form; see the module docstring."""

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

        for key, coeff in elem.items():
            ret += self.from_dict({_canonical_rc(key): coeff})
        return ret

    def mul(self, a, b):
        return self.from_dict(self._snap_qy(super().mul(a, b)))

    def rmul(self, a, b):
        return self.from_dict(self._snap_qy(super().rmul(a, b)))

    def from_dict(self, dct):
        elem = super().from_dict(dct)
        elem.ring = self
        return elem

from functools import cache

from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


class SchubertRCGraphRingElement(RCGraphRingElement):
    def to_free_algebra_element(self, basis=None):
        # from schubmult.free_algebra

        # result = sum([coeff * FSlideDual(*key.forest_weight) for key, coeff in self.items()])
        # if basis is None:
        #     return result
        # return result.change_basis(basis)
        pass


class SchubertRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = SchubertRCGraphRing._id
        SchubertRCGraphRing._id += 1
        self.dtype = type("SchubertRCGraphRingElement", (SchubertRCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbeasfrrtaffffystoa", self._ID))

    def new(self, x):
        return self.from_dict({x: 1}, snap=True)

    def __call__(self, x):
        return self.new(x)

    @property
    def zero_monom(self):
        return RCGraph([])

    def _snap(self, elem):
        # ret = self.zero
        # for key, coeff in elem.items():
        #     qy_code = key.snap_qy().length_vector
        #     the_rc = next(iter([rcc for rcc in RCGraph.all_rc_graphs(uncode(qy_code), len(key), weight=key.length_vector) if rcc.snap_qy().length_vector == qy_code]))
        #     # if the_rc is not None:
        #     ret += self.from_dict({the_rc: coeff})
        # return ret
        return elem

    @cache
    def _factor_schub(self, perm, length):
        from .bounded_rc_factor_algebra import BoundedRCFactorAlgebra
        r = BoundedRCFactorAlgebra()
        return r.schub_elem(perm, length)

    @cache
    def _factor_rc(self, rc):
        from .bounded_rc_factor_algebra import BoundedRCFactorAlgebra
        r = BoundedRCFactorAlgebra()
        return sum([r.from_dict({key: val}) for key, val in self._factor_schub(rc.perm, len(rc)).items() if r.key_to_rc_graph(key) == rc])

    def mul(self, a, b):
        result = 0
        for rc1, coeff1 in a.items():
            for rc2, coeff2 in b.items():
                rc = (self._factor_rc(rc1)*self._factor_rc(rc2)).to_rc_graph_ring_element()
                result += coeff1 * coeff2 * self._snap(rc)
        return result


    @cache
    def schubert_poly(self, perm, length):
        return self.from_dict(self.schub(perm, length))

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b), snap=True)

    def from_dict(self, dct, snap=False):
        elem = super().from_dict(dct)
        # elem.ring = self
        if snap:
            elem = self._snap(elem)
        # elem.ring = self
        return self.dtype({k: v for k, v in elem.items() if v != 0})


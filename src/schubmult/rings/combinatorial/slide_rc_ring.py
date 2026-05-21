from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


class SlideRCGraphRingElement(RCGraphRingElement):
    def to_free_algebra_element(self, basis=None):
        # from schubmult.free_algebra

        # result = sum([coeff * FSlideDual(*key.forest_weight) for key, coeff in self.items()])
        # if basis is None:
        #     return result
        # return result.change_basis(basis)
        pass


class SlideRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = SlideRCGraphRing._id
        SlideRCGraphRing._id += 1
        self.dtype = type("SlideRCGraphRingElement", (SlideRCGraphRingElement,), {"ring": self})

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
        ret = self.zero
        for key, coeff in elem.items():
            qy_code = key.snap_qy().length_vector
            the_rc = next(iter([rcc for rcc in RCGraph.all_rc_graphs(uncode(qy_code), len(key), weight=key.length_vector) if rcc.snap_qy().length_vector == qy_code]))
            # if the_rc is not None:
            ret += self.from_dict({the_rc: coeff})
        return ret

    @cache
    def _factor_schub(self, comp):
        from .bounded_rc_factor_algebra import BoundedRCFactorAlgebra
        r = BoundedRCFactorAlgebra()
        comp = tuple(comp)
        return r.from_dict({key: val for key, val in r.schub_elem(uncode(comp), len(uncode(comp))).items() if r.key_to_rc_graph(key).snap_qy().length_vector == comp})

    @cache
    def _factor_rc(self, rc):
        from .bounded_rc_factor_algebra import BoundedRCFactorAlgebra
        r = BoundedRCFactorAlgebra()
        #return r.from_dict({key: val for key, val in self._factor_schub(rc.snap_qy().length_vector).items() if r.key_to_rc_graph(key).length_vector == rc.length_vector})
        return r.from_rc_graph(rc, len(rc.perm))

    def mul(self, a, b):
        result = 0
        for rc1, coeff1 in a.items():
            for rc2, coeff2 in b.items():
                rc = (self._factor_rc(rc2)*self._factor_rc(rc1)).to_rc_graph_ring_element()
                result += coeff1 * coeff2 * self._snap(rc)
        return result


    @cache
    def slide_poly(self, comp):
        comp = tuple(comp)
        return self.from_dict({rc: 1 for rc in RCGraph.all_rc_graphs(uncode(comp), len(comp)) if rc.snap_qy().length_vector == comp})

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b), snap=True)

    def from_dict(self, dct, snap=False):
        elem = super().from_dict(dct)
        # elem.ring = self
        if snap:
            elem = self._snap(elem)
        # elem.ring = self
        return self.dtype({k: v for k, v in elem.items() if v != 0})


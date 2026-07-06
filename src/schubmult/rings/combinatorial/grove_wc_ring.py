from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.combinatorics.wc_graph import WCGraph

from .wc_graph_ring import WCGraphRing, WCGraphRingElement


def _canonical_rc(rc):
    # if rc.forest_weight == rc.length_vector:
    #     return rc
    # print(omega_park(rc.perm_word))
    # for rc0 in WCGraph.all_forest_rcs(rc.forest_weight):
    #     if rc0.forest_invariant == rc.forest_invariant and omega_park(tuple(reversed(rc0.perm_word))) == set(range(1, rc0.inv + 1)):
    #         return rc0
    #     print(rc0, rc0.forest_invariant, omega_park(tuple(reversed(rc0.perm_word))))
    # raise ValueError(f"Failed to find canonical RC graph for {rc}, which has forest weight {rc.forest_weight} and forest invariant {rc.forest_invariant}")
    st = [rc0 for rc0 in WCGraph.all_wc_graphs(uncode(rc.forest_weight), len(rc)) if rc0.omega_invariant[1] == rc.omega_invariant[1] and rc0.length_vector == rc.length_vector and rc0.forest_weight == rc.forest_weight and uncode(rc0.forest_weight) == rc0.perm]
    if len(st) != 1:
        raise ValueError(f"Failed to find unique canonical RC graph for {rc}, which has forest weight {rc.forest_weight} and omega invariant {rc.omega_invariant}. Candidates were: {st}")
    return st[0].resize(len(rc))
    # qy_rc = crc(rc)
    # the_real_qy = WCGraph.principal_rc(uncode(qy_rc.length_vector), len(rc))
    # return next(iter([rcc for rcc in WCGraph.all_wc_graphs(the_real_qy.perm, len(rc)) if crc(rcc) == the_real_qy and rcc.length_vector == rc.length_vector]))


class GroveWCGraphRingElement(WCGraphRingElement):

    def quasi_shift(self, i):
        result = self.ring.zero
        for key, coeff in self.items():
            if len(key) < i:
                result += coeff * self.ring(key)
                continue
            if key.length_vector[i - 1] != 0:
                continue
            result += coeff * self.ring(key.pull_out_row(i))
        return result

    def forest_trim(self, i):
        return self.divdiff(i).quasi_shift(i)

class DualGroveWCGraphRingElement(GroveWCGraphRingElement):
    pass

class GroveWCGraphRing(WCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = GroveWCGraphRing._id
        GroveWCGraphRing._id += 1
        self.dtype = type("GroveWCGraphRingElement", (GroveWCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbeasfrrtaffffsfsfystoa", self._ID))

    def new(self, x):
        return self.from_dict({x: 1}, snap=True)

    def __call__(self, x):
        return self.new(x)

    @property
    def zero_monom(self):
        return WCGraph([])

    def _snap(self, elem):
        ret = self.zero
        for key, coeff in elem.items():
            the_rc = _canonical_rc(key)
            # if the_rc is not None:
            ret += self.from_dict({the_rc: coeff})
        return ret

    def mul(self, a, b):
        result = self.zero
        for rc1, coeff1 in a.items():
            for rc2, coeff2 in b.items():
                if len(rc1) != len(rc2):
                    continue
                #     raise ValueError(f"Cannot multiply RC graphs of different lengths: {rc1} has length {len(rc1)}, while {rc2} has length {len(rc2)}")
                #     continue
                length = len(rc1)
                if length == 0:
                    # Identity times identity stays at the same ambient length (0).
                    # Avoid round-tripping through BoundedWCFactorAlgebra, which would
                    # otherwise inflate the empty RC graph to a nonzero ambient size.
                    result += coeff1 * coeff2 * self(rc1)
                    continue
                prd = (self._grove_bwc_lookup(rc1) * self._grove_bwc_lookup(rc2))
                result += coeff1 * coeff2 * self._snap(prd.to_wc_graph_ring_element().resize(length))
        return result

    @staticmethod
    @cache
    def _grove_bwc_lookup(rc):
        from schubmult.rings.combinatorial.bounded_wc_factor_algebra import BoundedWCFactorAlgebra
        r = BoundedWCFactorAlgebra()
        return r.from_wc_graph(rc, size=len(rc.perm))

    def coproduct_on_basis(self, a):
        addup = (self@self).zero
        for i in range(len(a) + 1):
            if i == 0:
                a1 = WCGraph([])
                a2 = a
            elif i == len(a):
                a1 = a
                a2 = WCGraph([])
            else:
                a1, a2 = a.vertical_cut(i)
            #if uncode(a1.forest_weight) == a1.perm or uncode(a2.forest_weight) == a2.perm:
            #if (self(a1) * self(a2))
            addup += self(a1) @ self(a2)
        return addup

    @cache
    def forest_poly(self, comp):
        comp = tuple(comp)
        return self.from_dict({rc: 1 for rc in WCGraph.all_wc_graphs(uncode(comp), len(comp)) if rc.forest_weight == comp})

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b), snap=True)

    def from_dict(self, dct, snap=False):
        elem = super().from_dict(dct)
        # elem.ring = self
        if snap:
            elem = self._snap(elem)
        # elem.ring = self
        return self.dtype({k: v for k, v in elem.items() if v != 0})


class DualGroveWCGraphRing(WCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = DualGroveWCGraphRing._id
        DualGroveWCGraphRing._id += 1
        self.dtype = type("DualGroveWCGraphRingElement", (DualGroveWCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbeasfrrtaffffsfsfystoa", self._ID))

    def new(self, x):
        return self.from_dict({x: 1}, snap=True)

    def __call__(self, x):
        return self.new(x)

    @property
    def zero_monom(self):
        return WCGraph([])

    def _snap(self, elem):
        ret = self.zero
        for key, coeff in elem.items():
            the_rc = _canonical_rc(key)
            if the_rc == key:
                ret += self.from_dict({the_rc: coeff})
        return ret

    def mul(self, a, b):
        return self.from_dict(super().mul(a, b), snap=True)


    def coproduct_on_basis(self, a):
        from schubmult.rings.free_algebra import ForestDual
        addup = (self@self).zero
        dr = GroveWCGraphRing()
        the_coprod = ForestDual(*a.forest_weight).coproduct()
        for (comp1, comp2), coeff in the_coprod.items():
            for rc1 in [rcc for rcc in WCGraph.all_wc_graphs(uncode(comp1), len(comp1)) if rcc.forest_weight == comp1]:
                for rc2 in [rcc for rcc in WCGraph.all_wc_graphs(uncode(comp2), len(comp2)) if rcc.forest_weight == comp2]:
                    if (dr(rc1) * dr(rc2)).almosteq(dr(a)):
                        addup += coeff * self(rc1) @ self(rc2)
        return addup

    @cache
    def forest_poly(self, comp):
        comp = tuple(comp)
        return self.from_dict({rc: 1 for rc in WCGraph.all_wc_graphs(uncode(comp), len(comp)) if rc.forest_weight == comp})

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b), snap=True)

    def from_dict(self, dct, snap=False):
        elem = super().from_dict(dct)
        # elem.ring = self
        if snap:
            elem = self._snap(elem)
        # elem.ring = self
        return self.dtype({k: v for k, v in elem.items() if v != 0})

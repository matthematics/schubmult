from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


def _canonical_rc(rc):
    # if rc.forest_weight == rc.length_vector:
    #     return rc
    # print(omega_park(rc.perm_word))
    # for rc0 in RCGraph.all_forest_rcs(rc.forest_weight):
    #     if rc0.forest_invariant == rc.forest_invariant and omega_park(tuple(reversed(rc0.perm_word))) == set(range(1, rc0.inv + 1)):
    #         return rc0
    #     print(rc0, rc0.forest_invariant, omega_park(tuple(reversed(rc0.perm_word))))
    # raise ValueError(f"Failed to find canonical RC graph for {rc}, which has forest weight {rc.forest_weight} and forest invariant {rc.forest_invariant}")
    st = [rc0 for rc0 in RCGraph.all_rc_graphs(uncode(rc.forest_weight), len(rc)) if rc0.omega_invariant[1] == rc.omega_invariant[1] and rc0.length_vector == rc.length_vector and rc0.forest_weight == rc.forest_weight]
    if len(st) != 1:
        raise ValueError(f"Failed to find unique canonical RC graph for {rc}, which has forest weight {rc.forest_weight} and omega invariant {rc.omega_invariant}. Candidates were: {st}")
    return st[0]
    # qy_rc = crc(rc)
    # the_real_qy = RCGraph.principal_rc(uncode(qy_rc.length_vector), len(rc))
    # return next(iter([rcc for rcc in RCGraph.all_rc_graphs(the_real_qy.perm, len(rc)) if crc(rcc) == the_real_qy and rcc.length_vector == rc.length_vector]))


class ForestRCGraphRingElement(RCGraphRingElement):
    def to_free_algebra_element(self, basis=None):
        from schubmult import ForestDual

        result = sum([coeff * ForestDual(*key.forest_weight) for key, coeff in self.items()])
        if basis is None:
            return result
        return result.change_basis(basis)


class ForestRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = ForestRCGraphRing._id
        ForestRCGraphRing._id += 1
        self.dtype = type("ForestRCGraphRingElement", (ForestRCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbeasfrrtaffffsfsfystoa", self._ID))

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
            the_rc = _canonical_rc(key)
            # if the_rc is not None:
            ret += self.from_dict({the_rc: coeff})
        return ret

    def mul(self, a, b):
        result = self.zero
        for rc1, coeff1 in a.items():
            for rc2, coeff2 in b.items():
                if len(rc1) != len(rc2):
                    raise ValueError(f"Cannot multiply RC graphs of different lengths: {rc1} has length {len(rc1)}, while {rc2} has length {len(rc2)}")
                    continue
                snap_size = max(len(rc1.perm.trimcode), len(rc2.perm.trimcode))
                rc1_resized = rc1.resize(snap_size)
                rc2_resized = rc2.resize(snap_size)
                prd = (self._forest_brc_lookup(rc1_resized) * self._forest_brc_lookup(rc2_resized))
                result += coeff1 * coeff2 * self._snap(prd.to_rc_graph_ring_element()).resize(len(rc1))
        return result


    @staticmethod
    @cache
    def _forest_brc(comp, length):
        from schubmult.rings.combinatorial.bounded_rc_factor_algebra import BoundedRCFactorAlgebra
        r = BoundedRCFactorAlgebra()
        #sb = r.full_schub_elem(uncode(comp), length)
        sb = r.schub_elem(uncode(comp), length)
        ret = r.from_tensor_dict({key: val for key, val in sb.items() if r.key_to_rc_graph(key).forest_weight == tuple(comp)}, size=length)
        if any(val < 0 for val in ret.values()):
            raise ValueError(f"Negative coefficient in forest BRC for {comp} of length {length}: {ret}")
        return ret

    @staticmethod
    @cache
    def _forest_brc_lookup(rc):
        brc = ForestRCGraphRing._forest_brc(rc.forest_weight, len(rc))
        return sum([coeff * brc.ring(key) for key, coeff in brc.items() if _canonical_rc(brc.ring.key_to_rc_graph(key)).resize(len(rc)) == _canonical_rc(rc)])

    def dual_product(self, a, b):
        pass

    @cache
    def forest_poly(self, comp):
        comp = tuple(comp)
        return self.from_dict({rc: 1 for rc in RCGraph.all_rc_graphs(uncode(comp), len(comp)) if rc.forest_weight == comp})

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b), snap=True)

    def from_dict(self, dct, snap=False):
        elem = super().from_dict(dct)
        # elem.ring = self
        if snap:
            elem = self._snap(elem)
        # elem.ring = self
        return self.dtype({k: v for k, v in elem.items() if v != 0})


from schubmult.combinatorics.indexed_forests import decreasing_labelings, omega_park, weak_composition_to_indfor, word_from_labeling
from schubmult.combinatorics.nilplactic import NilPlactic
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
    return next(iter([rc0 for rc0 in RCGraph.all_forest_rcs(rc.forest_weight) if rc0.forest_weight == rc0.length_vector and rc0.forest_invariant == rc.forest_invariant]))

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
            #if the_rc is not None:
            ret += self.from_dict({the_rc: coeff})
        return ret

    def mul(self, a, b):
        return self.from_dict(super().mul(a, b), snap=True)

    def rmul(self, a, b):
        return self.from_dict(super().rmul(a, b), snap=True)

    def from_dict(self, dct, snap=False):
        elem = super().from_dict(dct)
        #elem.ring = self
        if snap:
            elem = self._snap(elem)
        #elem.ring = self
        return self.dtype({k: v for k, v in elem.items() if v != 0})

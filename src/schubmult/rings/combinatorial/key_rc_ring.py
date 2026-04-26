from schubmult.combinatorics.rc_graph import RCGraph

from .rc_graph_ring import RCGraphRing, RCGraphRingElement


def _canonical_rc(rc):
    from schubmult import NilPlactic
    # not needed
    # if rc.extremal_weight != rc.perm.pad_code(len(rc)):
    #     return None
    # return rc#.to_highest_weight()[0]
    nilp, plac = NilPlactic.ed_column_insert_rsk(rc.perm_word, rc.compatible_sequence)
    for rc2 in RCGraph.all_key_rcs(rc.extremal_weight, weight=rc.length_vector):
        nilp2, plac2 = NilPlactic.ed_column_insert_rsk(rc2.perm_word, rc2.compatible_sequence)
        if plac == plac2:
            return rc
    raise ValueError("ca")
    # return None
    #return rc
    #the_set = {rc2 for rc2 in RCGraph.all_key_rcs(rc.extremal_weight, weight=rc.length_vector) if rc2.weight_tableau == rc.weight_tableau and rc.length_vector == rc2.length_vector}
    # # if len(the_set) > 1:
    # #     raise ValueError(f"Multiple RC graphs with same extremal weight and weight tableau found for {rc}, got {the_set}")
    #return next(iter(the_set))

class KeyRCGraphRingElement(RCGraphRingElement):
    def to_free_algebra_element(self, basis=None):
        from schubmult.rings.free_algebra import FreeAlgebra, KeyBasis
        KeyDual = FreeAlgebra(KeyBasis)
        result = sum([coeff * KeyDual(*key.extremal_weight) for key, coeff in self.items()])
        if basis is None:
            return result
        return result.change_basis(basis)

class KeyRCGraphRing(RCGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = KeyRCGraphRing._id
        KeyRCGraphRing._id += 1
        self.dtype = type("KeyRCGraphRingElement", (KeyRCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dasinkbeasfrrtaffffsfsfystoa", self._ID))

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
            if the_rc is not None:
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

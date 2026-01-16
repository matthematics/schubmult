from __future__ import annotations

from typing import TYPE_CHECKING

from schubmult.rings.schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement
from schubmult.schub_lib.bpd import BPD
from schubmult.schub_lib.perm_lib import Permutation

if TYPE_CHECKING:
    from schubmult.rings.rc_graph_ring import RCGraphRing, RCGraphRingElement


class BPDRingElement(SchubertMonomialRingElement):
    def to_rc_graph_ring_element(self, rc_ring: RCGraphRing | None = None) -> RCGraphRingElement:
        from schubmult.rings.rc_graph_ring import RCGraphRing

        if rc_ring is None:
            rc_ring = RCGraphRing()
        result_dict = {}
        for bpd, coeff in self.items():
            rc = bpd.to_rc_graph()
            result_dict[rc] = result_dict.get(rc, 0) + coeff
        return rc_ring.from_dict(result_dict)


class BPDRing(SchubertMonomialRing):
    _id = 0

    def __init__(self, genset=None, coeff_genset=None, *_, **__):
        # BPDs don't really use gensets since they're explicit graphical objects
        # But we need to initialize parent for BaseSchubertRing
        from schubmult.rings.variables import GeneratingSet
        from schubmult.symbolic import EXRAW

        if genset is None:
            genset = GeneratingSet("x")
        if coeff_genset is None:
            coeff_genset = GeneratingSet(None)

        # Initialize parent attributes directly instead of calling super().__init__
        # to avoid conflict with zero_monom property
        self._genset = genset
        self._coeff_genset = coeff_genset
        self.symbols = list(genset)
        self.domain = EXRAW
        self.dom = self.domain
        # Don't set self.zero_monom here - it's a property

        self._ID = BPDRing._id
        BPDRing._id += 1
        self.dtype = type("BPDRingElement", (BPDRingElement,), {"ring": self})


    @property
    def one(self) -> BPDRingElement:
        # Define the "one" element for BPDRing
        identity_graph = BPD.rothe_bpd(Permutation([]), 0)

        return self.from_dict({identity_graph: 1})

    @property
    def zero_monom(self):
        return BPD.rothe_bpd(Permutation([]), 0)

    def from_rc_graph_ring_element(self, rc_element: RCGraphRingElement) -> BPDRingElement:
        result_dict = {}
        for rc, coeff in rc_element.items():
            bpd = BPD.from_rc_graph(rc)
            result_dict[bpd] = result_dict.get(bpd, 0) + coeff
        return self.from_dict(result_dict)

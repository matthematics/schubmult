"""`ChuteMoveRing`: a `SchubertMonomialRing` whose basis elements are `ChuteMoveElement`s
(RC graphs marked with a set of simultaneous chute-move rows).
"""

from schubmult.combinatorics.chute_move_element import ChuteMoveElement
from schubmult.rings.combinatorial.schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class ChuteMoveRingElement(SchubertMonomialRingElement):
    """ChuteMoveRing elements are linear combinations of ChuteMoveElement basis elements."""

    # ----------------------
    # Presentation helpers
    # ----------------------



class ChuteMoveRing(SchubertMonomialRing):
    """The ring of `ChuteMoveElement`s; products use `ChuteMoveElement.product`."""

    _id = 0

    def __init__(self, *_, **__):
        self._ID = ChuteMoveRing._id
        ChuteMoveRing._id += 1
        self.dtype = type("ChuteMoveRingElement", (ChuteMoveRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkbesdgrrtystoa", self._ID))

    @property
    def zero_monom(self):
        return ChuteMoveElement([],[])

    @property
    def one(self):
        # Define the "one" element for ChuteMoveRing
        identity_graph = ChuteMoveElement([],[])
        return self.from_dict({identity_graph: 1})

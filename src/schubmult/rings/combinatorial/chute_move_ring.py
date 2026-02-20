from schubmult.combinatorics.chute_move_element import ChuteMoveElement
from schubmult.rings.combinatorial.schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class ChuteMoveRingElement(SchubertMonomialRingElement):
    """
    ChuteMoveRing elements are linear combinations of ChuteMoveElement basis elements.

    The product % is the polynomial product. Currently only defined when the right side
    is a dominant RC graph.

    The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
    to do this directly.

    The product * is well defined for any pair of RC graphs and is the dual product.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------



class ChuteMoveRing(SchubertMonomialRing):
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

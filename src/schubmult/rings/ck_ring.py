from functools import cache

from sympy import pretty_print

from schubmult import ASx
from schubmult.perm_lib import uncode
from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.crystal_graph import CrystalGraphTensor
from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.rings.plactic import NilPlactic, Plactic
from schubmult.symbolic import S, sympy_Mul

from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from .rc_graph import RCGraph
from .schubert_ring import Sx

#weight wt
# yw highest weight
# u # yv
# yv highest weight


class CoxeterKnuthPrintingTerm(TypedPrintingTerm):
    pass

class CoxeterKnuthRingElement(CrystalGraphRingElement):
    """
    CoxeterKnuthRing elements are linear combinations of NilPlactic basis elements.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------
    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [
            self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k))
            for k in self.keys()
        ]

    def __eq__(self, other):
        return type(self) is type(other) and dict(self) == dict(other)

class CoxeterKnuthKey(Plactic):

    def __init__(self, p_tableau, weight_tableau, length):
        self._p_tableau = p_tableau
        self._word = weight_tableau._word
        self._weight_tableau = weight_tableau
        self.length = length

    def crystal_length(self) -> int:
        return self.length

    def _sympystr(self, printer):
        return printer._print(CrystalGraphTensor(self._p_tableau, self._weight_tableau))

    def _pretty(self, printer):
        return printer._print(CrystalGraphTensor(self._p_tableau, self._weight_tableau))

    @property
    def rc_graph(self):
        rc = self._p_tableau.hw_rc(self.length)
        _, raise_seq = self._weight_tableau.to_highest_weight(length=self.length)
        rc = rc.reverse_raise_seq(raise_seq)
        return rc
    
    @classmethod
    def from_rc_graph(cls, rc: RCGraph):
        return cls(rc.p_tableau, rc.weight_tableau, len(rc))

    def monk_insert(self, simple_box, monk_top, retries = 6, target_perms=None):

        """
        Candidate associative 'Monk' insertion for this CoxeterKnuthKey.

        Parameters
        - simple_box: either a Plactic instance representing a one-box tableau,
                      or an int giving the entry of that single box.

        Behavior (candidate):
        - Insert the simple_box into the weight tableau using rs_insert (unique).
        - Attempt to update the P-tableau by applying the corresponding simple
          action to the underlying RC-graph. We try several RCGraph APIs
          (act, iterative_act, prod_with_rc) in order and collect resulting
          RC-graphs. For each resulting RC-graph we produce a new
          CoxeterKnuthKey whose weight tableau is the rs-inserted tableau and
          whose P-tableau is taken from the RC-graph's p_tableau.
        - Return: list of CoxeterKnuthKey instances (may be empty if no
          candidate RC-graphs were produced).

        Note: this is a conservative "candidate" implementation; it does not
        attempt to prove Monk-associativity. It simply follows the natural
        route: shift/update weight by rs_insert and propagate the action to
        the rc_graph, pulling back the P-tableau from any produced RC-graphs.
        """
        # Normalize simple_box to a single integer entry
        assert simple_box <= monk_top, "Cannot insert box larger than monk_top"
        box_val = simple_box
        
        # 2) attempt to produce candidate updated P-tableaux by acting on rc_graph
        rc0 = self.rc_graph

        if len(rc0) < monk_top:
            rc0 = rc0.resize(monk_top)

        poly = Sx(~self._p_tableau.perm)
        poly *= Sx(uncode(([0]* (monk_top - 1)) + [1]))
        if target_perms is not None:
            target_perms = target_perms.intersection(set(poly.keys()))
        else:
            target_perms = set(poly.keys())
        # try rc.act(box_val)
        candidates = set()
        weight = [*rc0.length_vector]
        weight[box_val - 1] += 1
        target_weight = tuple(weight)
        for perm in target_perms:
            candidates.update(RCGraph.all_rc_graphs(perm, len(rc0), weight=target_weight))
        weight2 = self._weight_tableau.rs_insert(box_val)
        results = set()
        for rc in candidates:
            if rc.weight_tableau == weight2:
                assert rc.perm.inv == self._p_tableau.perm.inv + 1
                results.add(rc)
        if len(results) > 0:
            return results
        if retries > 0:
            target_perms2 = set()
            for perm in target_perms:
                target_perms2.add(uncode([0, *perm.trimcode]))
            res =  self.shiftup(1).monk_insert(box_val + 1, max(box_val + 1, len((~self._p_tableau.perm).trimcode) + 1), retries - 1, target_perms=target_perms2)
            if res is not None:
                for rc in res:
                    assert rc.perm.inv == self._p_tableau.perm.inv + 1
                    results.add(rc.rowrange(1))
                if len(results) > 0:
                    return results
        return None
        # try multiplying by a one-row RCGraph corresponding to the box (fallback)
        #return CoxeterKnuthKey.from_rc_graph(RCGraph([(),*self.rc_graph.shiftup(1)])).monk_insert(box_val + 1).shiftup(-1).normalize()
    
    def shiftup(self, k1):
        return CoxeterKnuthKey(NilPlactic.from_word(self._p_tableau.shiftup(k1)), self._weight_tableau.shiftup(k1), self.length + k1)



class CoxeterKnuthRing(CrystalGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = CoxeterKnuthRing._id
        CoxeterKnuthRing._id += 1
        self.dtype = type("CoxeterKnuthRingElement", (CoxeterKnuthRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtystoa", "poing", self._ID))

    @property
    def zero_monom(self):
        return (NilPlactic(), 0)


    def printing_term(self, key):
        return CoxeterKnuthPrintingTerm(key)

    # def dtype(self):
    #     elem = CoxeterKnuthRingElement()
    #     elem.ring = self
    #     return elem

    def from_dict(self, dct):
        elem = self.dtype()
        elem.update(dct)
        return elem

    def __call__(self, key):
        return self.from_dict({CoxeterKnuthKey(*key): 1})

    def mul(self, a, b):
        # a, b are CoxeterKnuthRingElemen
        if isinstance(b, CoxeterKnuthRingElement):
            result_dict = {}
            for g1, c1 in a.items():
                for g2, c2 in b.items():
                    # CoxeterKnuth.prod_with_rc returns a dict {CoxeterKnuth: coeff}
                    prod = g1.rc_graph.prod_with_rc(g2.rc_graph)
                    for g3, c3 in prod.items():
                        result_dict[CoxeterKnuthKey(g3.p_tableau,g3.weight_tableau,len(g3))] = result_dict.get(CoxeterKnuthKey(g3.p_tableau,g3.weight_tableau,len(g3)), 0) + c1 * c2 * c3
            # result_dict = {k: v * b for k, v in a.items()}
        return self.from_dict(result_dict)

    def __eq__(self, other):
        return type(self) is type(other) and self._ID == other._ID

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        # Define the "one" element for CoxeterKnuthRing
        identity_graph = (NilPlactic(), Plactic(), 0)
        return self.from_dict({CoxeterKnuthKey(*identity_graph): 1})


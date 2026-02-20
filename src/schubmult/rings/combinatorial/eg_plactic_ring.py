from sympy import Tuple

from schubmult.combinatorial_reps.nilplactic import NilPlactic

# from schubmult.combinatorial_reps.nilplactic import NilPlactic
# from schubmult.combinatorial_reps.plactic import Plactic
from schubmult.combinatorial_reps.rc_graph import RCGraph
from schubmult.symbolic import S

from ..printing import PrintingTerm, TypedPrintingTerm
from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from .rc_graph_ring import RCGraphRing

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class EGPlacticPrintingTerm(PrintingTerm):
    is_commutative = False
    precedence = 50

    def __new__(cls, k):
        return EGPlacticPrintingTerm.__xnew_cached__(cls, k)

    @staticmethod
    def __xnew__(_class, k):
        obj = PrintingTerm.__new__(_class, k, None, None)
        obj._key = k
        return obj

    def __hash__(self):
        return hash((self._key, "barclic"))

    @staticmethod
    def __xnew_cached__(_class, k):
        return EGPlacticPrintingTerm.__xnew__(_class, k)

    def _sympystr(self, printer):
        return " # ".join([printer._print(self._key[i]) for i in range(len(self._key))])

    def _pretty(self, printer):
        return printer._print_TensorProduct(Tuple(TypedPrintingTerm(self._key[0]), TypedPrintingTerm(self._key[1])))

    def _latex(self, printer):
        return printer._print_TensorProduct(Tuple(TypedPrintingTerm(self._key[0]), TypedPrintingTerm(self._key[1])))


class EGPlacticRingElement(CrystalGraphRingElement):
    """
    EGPlacticRing elements are linear combinations of RCGraph basis elements.

    The product % is the polynomial product. Currently only defined when the right side
    is a dominant RC graph.

    The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
    to do this directly.

    The product * is well defined for any pair of RC graphs and is the dual product.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------

    def to_rc_graph_ring_element(self):
        r = RCGraphRing()
        return r.from_dict({self.ring.key_to_rc_graph(k): v for k, v in self.items()})

    def __mod__(self, other):
        """
        Polynomial product: self % other.
        Currently only defined when `other` is a dominant RC graph.
        """
        return self.ring.rc_product(self, other)

    def raising_operator(self, index):
        """
        Linear extension of RCGraph.raising_operator:
        Apply raising_operator(index) to every basis RCGraph in self, collect results.
        Returns an EGPlacticRingElement (possibly zero).
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc = rc_graph.raising_operator(index)
            if new_rc is not None:
                res += coeff * self.ring(new_rc)
        return res

    def lowering_operator(self, index):
        """
        Linear extension of RCGraph.lowering_operator.
        """
        res = self.ring.zero
        for rc_graph, coeff in self.items():
            new_rc = rc_graph.lowering_operator(index)
            if new_rc is not None:
                res += coeff * self.ring(new_rc)
        return res

    def phi(self, index):
        """
        phi(element) := max_{basis rc in supp(element)} phi(rc)
        If element is zero, returns 0.
        """
        if len(self) == 0:
            return 0
        m = 0
        for rc_graph in self.keys():
            try:
                v = rc_graph.phi(index)
            except Exception:
                v = 0
            if v > m:
                m = v
        return m

    def epsilon(self, index):
        """
        epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)
        """
        if len(self) == 0:
            return 0
        m = 0
        for rc_graph in self.keys():
            try:
                v = rc_graph.epsilon(index)
            except Exception:
                v = 0
            if v > m:
                m = v
        return m

    def crystal_length(self):
        """
        Use maximum crystal length of basis graphs in support (0 for the zero element).
        """
        if len(self) == 0:
            return 0
        return max(getattr(rc_graph, "crystal_length", lambda: len(rc_graph))() for rc_graph in self.keys())

    def almosteq(self, other):
        return all(v == 0 for v in (self - other).values())

    def to_highest_weight(self):
        """
        Iteratively raise the element until no further raising is possible.
        Returns (highest_weight_element, raise_seq).

        Behavior notes:
        - This is the natural linear-extension of CrystalGraph.to_highest_weight.
        - The returned `highest_weight_element` is an EGPlacticRingElement.
        - raise_seq is the sequence of row indices applied (in order).
        """
        rc_elem = self.ring.zero
        for rc, coeff in self.items():
            rc_elem += coeff * self.ring(rc.to_highest_weight()[0])
        return rc_elem, None

    def to_lowest_weight(self):
        """
        Iteratively raise the element until no further raising is possible.
        Returns (highest_weight_element, raise_seq).

        Behavior notes:
        - This is the natural linear-extension of CrystalGraph.to_highest_weight.
        - The returned `highest_weight_element` is an EGPlacticRingElement.
        - raise_seq is the sequence of row indices applied (in order).
        """
        rc_elem = self.ring.zero
        for rc, coeff in self.items():
            rc_elem += coeff * self.ring(rc.to_lowest_weight()[0])
        return rc_elem, None

    def reverse_raise_seq(self, raise_seq):
        """
        Apply lowering_operator in reverse order to `raise_seq`.
        If the path dies (result is zero), return None (mirrors scalar behavior).
        """
        rc_elem = self
        for row in reversed(raise_seq):
            rc_elem = rc_elem.lowering_operator(row)
            if rc_elem is None or len(rc_elem) == 0:
                return None
        return rc_elem


class EGPlacticRing(CrystalGraphRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = EGPlacticRing._id
        EGPlacticRing._id += 1
        self.dtype = type("EGPlacticRingElement", (EGPlacticRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dsdfginkbeafrrtystoa", self._ID))

    def from_rc_graph(self, rc):
        nilp, plac = rc.hw_tab_rep()
        nilp = NilPlactic.from_word(rc.perm_word)
        return self.from_dict({((nilp, len(rc)), plac): S.One})

    @property
    def zero_monom(self):
        return RCGraph([])

    def from_rc_graph_ring_element(self, elem):
        return self.from_dict({((NilPlactic.from_word(rc.perm_word), len(rc)), rc.hw_rc_rep()[1]): coeff for rc, coeff in elem.items() if coeff != 0})

    def from_dict(self, d):
        elem = EGPlacticRingElement()
        elem.update(d)
        elem.ring = self
        return elem

    def printing_term(self, key):
        return EGPlacticPrintingTerm(key)

    @property
    def one(self):
        # Define the "one" element for EGPlacticRing
        identity_graph = RCGraph()
        return self.from_rc_graph(identity_graph)

    def mul_pair(self, key1, key2, check=False):  # noqa: ARG002
        amt_to_bump = max(key2[0][1], len(key2[0][0].perm))
        # r = RCGraphRing()
        nilp_set = {key1[0][0]}
        self_len = key1[0][1]
        for a in range(amt_to_bump):
            new_nilp_set = set()
            for nilp in nilp_set:
                st = nilp.right_zero_act(self_len + a)
                new_nilp_set.update(st)
            nilp_set = new_nilp_set

        # rc_elem = r.zero
        ret = self.zero
        for nilp in nilp_set:
            # new_nilp = NilPlactic.from_word([*nilp.column_word, *[a + self_len for a in key2[0][0].column_word]])
            # negate?
            hw_plac1, raise_seq1 = key1[1].to_highest_weight(length=key1[0][1])
            hw_plac2, raise_seq2 = key2[1].to_highest_weight(length=key2[0][1])
            new_nilp, plac_elem_hw = NilPlactic.ed_column_insert_rsk(
                [*nilp.column_word, *[a + self_len for a in key2[0][0].column_word]],
                sorted([*hw_plac1.row_word] + [a + self_len for a in hw_plac2.row_word]),
            )
            if new_nilp.perm.inv == len(new_nilp.column_word) and len(new_nilp.perm.trimcode) <= self_len + key2[0][1]:
                plac_elem = plac_elem_hw.reverse_raise_seq(raise_seq1)
                plac_elem = plac_elem.reverse_raise_seq([a + self_len for a in raise_seq2])
                ret += self((new_nilp, self_len + key2[0][1]), plac_elem)
        return ret

    def key_to_rc_graph(self, key):
        _, raise_seq = key[1].to_highest_weight(length=key[0][1])
        return key[0][0].hw_rc(key[0][1]).reverse_raise_seq(raise_seq)

    def __call__(self, rc, tab):
        # if not rc.is_highest_weight:
        #     raise ValueError("rc must be highest weight")
        return self.from_dict({(rc, tab): S.One})

    def mul(self, a, b):
        result = self.zero
        for key1, coeff1 in a.items():
            for key2, coeff2 in b.items():
                result += coeff1 * coeff2 * self.mul_pair(key1, key2)
        return result

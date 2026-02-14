from sympy import Tuple

from schubmult.schub_lib.plactic import Plactic
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.symbolic import S
from schubmult.utils.perm_utils import little_bump

from .abstract_schub_poly import AbstractSchubPoly, TypedPrintingTerm
from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from .rc_graph_ring import RCGraphRing

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class EGPlacticPrintingTerm(AbstractSchubPoly):
    is_commutative = False
    precedence = 50

    def __new__(cls, k):
        return EGPlacticPrintingTerm.__xnew_cached__(cls, k)

    @staticmethod
    def __xnew__(_class, k):
        obj = AbstractSchubPoly.__new__(_class, k, None, None)
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
        return hash(("Dsdfginkberrtystoa", self._ID))

    def from_rc_graph(self, rc):
        return self.from_dict({rc.hw_tab_rep(): S.One})

    @property
    def zero_monom(self):
        return RCGraph([])

    def from_rc_graph_ring_element(self, elem):
        return self.from_dict({rc.hw_tab_rep(): coeff for rc, coeff in elem.items() if coeff != 0})

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

    @staticmethod
    def little_zero(word, length):
        from schubmult import Permutation
        perm = Permutation.ref_product(*word)
        if len(perm.trimcode) < length:
            return tuple(word)
        if len(perm.trimcode) > length:
            raise ValueError("Word is too long for the specified length")
        new_word = [*word]
        while len(perm.trimcode) >= length:
            d = len(perm.trimcode)
            old_word = new_word
            new_word = little_bump(new_word, d, d + 1)
            if new_word == old_word:
                raise ValueError(f"Word cannot be bumped further {word=} {new_word} {length=}")
            perm = Permutation.ref_product(*new_word)
        return tuple(new_word)


    @staticmethod
    def  _right_zero_act(word, length):
        from schubmult import ASx, Permutation, uncode
        perm = Permutation.ref_product(*word)
        up_perms = ASx(perm, length) * ASx(uncode([0]), 1)

        word_set = set()
        word = tuple(word)
        for perm1, _ in up_perms.keys():
            for rc in RCGraph.all_hw_rcs(perm1, length + 1):
                new_word = EGPlacticRing.little_zero(rc.perm_word, length + 1)
                if new_word == word:
                    word_set.add(rc.perm_word)
        return word_set

    def mul_pair(self, key1, key2, check=True):
        from schubmult import RCGraphRing

        amt_to_bump = max(len(key2[0]), len(key2[0].perm))
        r = RCGraphRing()
        rc_elem = r(key1[0])
        self_len = len(key1[0])
        for a in range(amt_to_bump):
            new_rc_elem = r.zero
            for rcc, coeff in rc_elem.items():
                word = rcc.perm_word
                st = EGPlacticRing._right_zero_act(word, len(rcc))
                for tup in st:
                    seq = []
                    last_spot = 0
                    last_elem = -1000
                    for a in tup:
                        if a > last_elem:
                            last_elem = a
                            last_spot += 1
                        seq.append(last_spot)
                        last_elem = a
                    new_rc_elem += coeff * r(RCGraph.from_reduced_compatible(tup, seq).resize(len(rcc) + 1))
            rc_elem = new_rc_elem

        old_rc_elem = rc_elem
        rc_elem = r.zero
        for rc, coeff in old_rc_elem.items():
            new_rc = RCGraph([*rc[:self_len],*key2[0].shiftup(self_len)]).resize(self_len + len(key2[0]))
            if new_rc.is_valid and len(new_rc.perm.trimcode) <= self_len + len(key2[0]):
                rc_elem += coeff * r(new_rc)
        if check:
            real_rc_elem = r(key1[0]) * r(key2[0])
            assert real_rc_elem.almosteq(rc_elem), f"Failed on {key1} and {key2} dingbat\n{real_rc_elem}\nstinkbat\n{rc_elem}\n{real_rc_elem - rc_elem}"
        ret = self.zero
        for new_rc, coeff in rc_elem.items():
            hw_rc, raise_seq3 = new_rc.to_highest_weight()
            _, raise_seq1 = key1[1].to_highest_weight(length=len(key1[0]))
            _, raise_seq2 = key2[1].to_highest_weight(length=len(key2[0]))
            plac_elem = Plactic.yamanouchi([a for a in hw_rc.length_vector if a != 0])
            plac_elem = plac_elem.reverse_raise_seq(raise_seq3)
            plac_elem = plac_elem.reverse_raise_seq(raise_seq1)
            plac_elem = plac_elem.reverse_raise_seq([a + len(key1[0]) for a in raise_seq2])
            ret += coeff * self(hw_rc, plac_elem)
        return ret

    def key_to_rc_graph(self, key):
        _, raise_seq = key[1].to_highest_weight(length=len(key[0]))
        return key[0].reverse_raise_seq(raise_seq)

    def __call__(self, rc, tab):
        if not rc.is_highest_weight:
            raise ValueError("rc must be highest weight")
        return self.from_dict({(rc, tab): S.One})

    def mul(self, a, b):
        result = self.zero
        for key1, coeff1 in a.items():
            for key2, coeff2 in b.items():
                result += coeff1 * coeff2 * self.mul_pair(key1, key2)
        return result

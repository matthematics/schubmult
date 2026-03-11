from schubmult.combinatorics.nilplactic import NilPlactic
from schubmult.symbolic import S

from ..base_ring import BaseRing, BaseRingElement
from ..printing import PrintingTerm


class EGPrintingTerm(PrintingTerm):
    is_commutative = False
    precedence = 50

    def __new__(cls, k):
        return EGPrintingTerm.__xnew_cached__(cls, k)

    @staticmethod
    def __xnew__(_class, k):
        obj = PrintingTerm.__new__(_class, k, None, None)
        obj._key = k
        return obj

    def __hash__(self):
        return hash((self._key, "barclic"))

    @staticmethod
    def __xnew_cached__(_class, k):
        return EGPrintingTerm.__xnew__(_class, k)

    def _sympystr(self, printer):
        return printer._print(self._key)

    def _pretty(self, printer):
        return printer._print(self._key)

    def _latex(self, printer):
        return printer._print(self._key)


class EGRingElement(BaseRingElement):
    """
    EGRing elements are linear combinations of RCGraph basis elements.

    The product % is the polynomial product. Currently only defined when the right side
    is a dominant RC graph.

    The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
    to do this directly.

    The product * is well defined for any pair of RC graphs and is the dual product.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------



class EGRing(BaseRing):
    _id = 0

    def __init__(self, *_, **__):
        self._ID = EGRing._id
        EGRing._id += 1
        self.dtype = type("EGRingElement", (EGRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dsdfginkberysfhafrrtystoa", self._ID))

    def from_rc_graph(self, rc):
        nilp, plac = rc.hw_tab_rep()
        nilp = NilPlactic.from_word(rc.perm_word)
        return self.from_dict({(nilp, len(rc)): S.One})

    @property
    def zero_monom(self):
        return (NilPlactic(), 0)

    def from_rc_graph_ring_element(self, elem):
        return self.from_dict({(NilPlactic.from_word(rc.perm_word), len(rc)): coeff for rc, coeff in elem.items() if coeff != 0})

    def from_dict(self, d):
        elem = EGRingElement()
        elem.update(d)
        elem.ring = self
        return elem

    def printing_term(self, key):
        return EGPrintingTerm(key)

    @property
    def one(self):
        # Define the "one" element for EGRing
        return {self.zero_monom: S.One}

    def mul_pair(self, key1, key2, check=False):  # noqa: ARG002
        amt_to_bump = max(key2[1], len(key2[0].perm))
        # r = RCGraphRing()
        nilp_set = {key1[0]}
        self_len = key1[1]
        for a in range(amt_to_bump):
            new_nilp_set = set()
            for nilp in nilp_set:
                st = nilp.right_zero_act(self_len + a)
                new_nilp_set.update(st)
            nilp_set = new_nilp_set

        ret = self.zero
        for nilp in nilp_set:
            new_nilp, plac_elem_hw = NilPlactic.ed_insert(
                *nilp.column_word, *[a + self_len for a in key2[0].column_word],
            )
            if new_nilp.perm.inv == len(new_nilp.column_word) and len(new_nilp.perm.trimcode) <= self_len + key2[0][1]:
                ret += self(new_nilp, self_len + key2[1])
        return ret

    def __call__(self, nilp, length):
        return self.from_dict({(nilp, length): S.One})

    def mul(self, a, b):
        result = self.zero
        for key1, coeff1 in a.items():
            for key2, coeff2 in b.items():
                result += coeff1 * coeff2 * self.mul_pair(key1, key2)
        return result

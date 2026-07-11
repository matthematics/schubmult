
from schubmult.combinatorics.wc_graph import WCGraph
from schubmult.rings.free_algebra import FreeAlgebra, GrothendieckBasis, WordBasis
from schubmult.symbolic import S, sympify_sympy, sympy_Mul

from ..schubert.grothendieck_ring import Gx

#from .crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from .schubert_monomial_ring import SchubertMonomialRing, SchubertMonomialRingElement

# from .crystal_graph_ring import CrystalTensorRing

# weight wt
# yw highest weight
# u # yv
# yv highest weight


class WCGraphRingElement(SchubertMonomialRingElement):
    """
    WCGraphRing elements are linear combinations of WCGraph basis elements.

    The product % is the polynomial product. Currently only defined when the right side
    is a dominant RC graph.

    The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
    to do this directly.

    The product * is well defined for any pair of RC graphs and is the dual product.
    """

    # ----------------------
    # Presentation helpers
    # ----------------------

    def as_terms(self):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        # Keep WCGraph([]) as a basis monomial rather than collapsing to scalar 1.
        return [sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in self.keys()]

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympify_sympy(S.Zero)]
        # Keep WCGraph([]) as a basis monomial rather than collapsing to scalar 1.
        return [sympy_Mul(sympify_sympy(self[k]), self.ring.printing_term(k)) for k in sorted(self.keys())]

    def almosteq(self, other):
        return all(v == 0 for v in (self - other).values())


    def zero_out_last_row(self):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            if len(wc_graph) == 0:
                res += coeff * self.ring(wc_graph)
            elif len(wc_graph[-1]) == 0:
                res += coeff * self.ring(wc_graph.zero_out_last_row())
        return res

    def resize(self, n):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            res += coeff * self.ring(wc_graph.resize(n))
        return res

    def clip(self, n):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            res += coeff * self.ring(wc_graph.vertical_cut(n)[0])
        return res

    def transpose(self, length):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            res += coeff * self.ring(wc_graph.transpose(length))
        return res

    def project(self):
        return self.ring.from_free_algebra_element(self.to_free_algebra_element())

    def trim_operator(self, i):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            res += coeff * self.ring.trim_operator(i, wc_graph)
        return res

    def double_elem_sym_squash(self, weight, yvars, zvars):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            res += coeff * wc_graph.double_elem_sym_squash(weight, yvars, zvars)
        return res

    def full_double_elem_sym_squash(self, p, yvars, zvars):
        res = self.ring.zero
        for wc_graph, coeff in self.items():
            res += coeff * wc_graph.full_double_elem_sym_squash(p, yvars, zvars)
        return res

    def grass_coaction(self):
        ret = self.ring.zero @ self.ring.zero
        for rc, coeff in self.items():
            ret += coeff * self.ring.grass_coaction(rc)
        return ret

    @property
    def broadcast(self):
        class BroadcastWrapper:
            def __init__(self, elem):
                self.elem = elem

            def __getattr__(self, attr):
                def method(*args, **kwargs):
                    res = self.elem.ring.zero
                    for wc_graph, coeff in self.elem.items():
                        func_result = getattr(wc_graph, attr)(*args, **kwargs)
                        res += coeff * self.elem.ring(func_result)
                    return res
                return method
        return BroadcastWrapper(self)


class WCGraphRing(SchubertMonomialRing):
    _id = 0

    def __call__(self, x):
        if isinstance(x, WCGraphRingElement):
            return self.from_dict(x)
        return self.new(x)

    def new(self, x):
        return self.from_dict({x: 1})

    def __init__(self, *_, **__):
        self._ID = WCGraphRing._id
        WCGraphRing._id += 1
        WCGraphRing._beta = Gx._beta
        self.dtype = type("WCGraphRingElement", (WCGraphRingElement,), {"ring": self})

    def __hash__(self):
        return hash(("Dinkberrtystoa", self._ID))

    def _weight_key(self, rc):
        """Return the canonical RC graph for *rc*'s (perm, length_vector) class."""
        return min(WCGraph.all_wc_graphs(rc.perm, len(rc), weight=rc.length_vector))

    def from_dict(self, dct):
        return super().from_dict(dct)

    def monomial(self, *tup):
        elem = self.one
        if len(tup) <= 1:
            for a in tup:
                elem = elem * self(WCGraph.one_row(a))
            return elem
        mid = len(tup) // 2
        return self.monomial(*tup[:mid]) * self.monomial(*tup[mid:])

    def elem_sym(self, descent, weight):
        from schubmult import uncode

        return self.from_dict(dict.fromkeys(WCGraph.all_wc_graphs(uncode([0] * (descent - sum(weight)) + [1] * sum(weight)), len(weight), weight=weight), 1))

    @property
    def zero_monom(self):
        return WCGraph([])

    def from_free_algebra_element(self, elem):
        wordelem = elem.change_basis(WordBasis)
        result = self.zero
        for word, coeff in wordelem.items():
            result += coeff * self.monomial(*word)
        return result

    def groth(self, perm, n=None):
        """
        Return the WCGraphRing element corresponding to the Schubert polynomial
        indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).
        """
        if n is None:
            n = len(perm.trimcode)
        return self.from_dict(dict.fromkeys(WCGraph.all_wc_graphs(perm, n), S.One))

    def to_free_algebra_element(self, elem, basis=GrothendieckBasis):
        dual_groth = FreeAlgebra(GrothendieckBasis)
        result = 0
        for k, v in elem.items():
            result += v * dual_groth((k.perm, len(k)))
        if basis != GrothendieckBasis:
            result = result.change_basis(basis)
        return result

    def trim_operator(self, i, rc):
        res = self(rc).divdiff(i)
        ret = self.zero
        for rc0, coeff in res.items():
            if len(rc0[i - 1]) != 0:
                continue
            ret += coeff * self(rc0.pull_out_row(i, keep_size=True))
        return ret



    @property
    def one(self):
        # Define the "one" element for WCGraphRing
        identity_graph = WCGraph()
        return self.from_dict({identity_graph: 1})

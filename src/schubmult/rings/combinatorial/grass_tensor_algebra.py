from __future__ import annotations

from sympy import Tuple

from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.base_ring import BaseRing, BaseRingElement
from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRing, RCGraphRingElement
from schubmult.rings.printing import PrintingTerm, TypedPrintingTerm
from schubmult.symbolic import S


def _is_full_grassmannian_rc(rc: RCGraph) -> bool:
    try:
        attr = rc.is_full_grassmannian
    except AttributeError:
        attr = None
    if isinstance(attr, bool):
        return attr
    return rc.perm.inv == 0 or rc.perm.descents() == {len(rc) - 1}

def _ensure_valid_key(key):
    if not isinstance(key, tuple):
        raise TypeError(f"Expected tuple key, got {type(key)}")
    for rc in key:
        if not isinstance(rc, RCGraph):
            raise TypeError(f"Key factors must be RCGraph, got {type(rc)}")
        if not _is_full_grassmannian_rc(rc):
            raise ValueError(f"Key factors must be full Grassmannian RC graphs, got {rc}")

def _descent_of_grass(rc: RCGraph) -> int:
    if rc.perm.inv == 0:
        return -1
    descs = rc.perm.descents()
    if len(descs) == 0:
        return -1
    return max(descs)


def _last_descent_size(rc: RCGraph) -> int:
    """Return max descent + 1, or 0 for identity."""
    if rc.perm.inv == 0:
        return 0
    descs = rc.perm.descents()
    if len(descs) == 0:
        return 0
    return max(descs) + 1


class GrassTensorPrintingTerm(PrintingTerm):
    is_commutative = False
    precedence = 50

    def __new__(cls, key):
        return GrassTensorPrintingTerm.__xnew_cached__(cls, key)

    @staticmethod
    def __xnew__(_class, key):
        obj = PrintingTerm.__new__(_class, key, None, None)
        obj._key = key
        return obj

    @staticmethod
    def __xnew_cached__(_class, key):
        return GrassTensorPrintingTerm.__xnew__(_class, key)

    def __hash__(self):
        return hash((self._key, "GrassTensorPrintingTerm"))

    def _sympystr(self, printer):
        if len(self._key) == 0:
            return printer._print(S.One)
        return " # ".join(printer._print(TypedPrintingTerm(factor)) for factor in self._key)

    def _pretty(self, printer):
        if len(self._key) == 0:
            return printer._print(S.One)
        return printer._print_TensorProduct(Tuple(*[TypedPrintingTerm(factor) for factor in self._key]))

    def _latex(self, printer):
        if len(self._key) == 0:
            return printer._print(S.One)
        return printer._print_TensorProduct(Tuple(*[TypedPrintingTerm(factor) for factor in self._key]))


class GrassTensorAlgebraElement(BaseRingElement):
    """Element of GrassTensorAlgebra: finite linear combinations of Grass tensors."""

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else self[k] * self.ring.printing_term(k) for k in self.keys()]

    def to_rc_graph_ring_element(self) -> RCGraphRingElement:
        r = RCGraphRing()
        result = r.zero
        for key, coeff in self.items():
            if coeff == 0:
                continue
            rc = self.ring.key_to_rc_graph(key)
            result += coeff * r(rc)
        return result

class GrassTensorAlgebra(BaseRing):
    """Tensor-like algebra on tuples of full Grassmannian RC graphs.

    Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
    RC graph. Simplification rules
    """

    _id = 0
    zero_monom = ()

    def __init__(self, domain=None):
        super().__init__(domain=domain)
        self._ID = GrassTensorAlgebra._id
        GrassTensorAlgebra._id += 1
        self.dtype = type("GrassTensorAlgebraElement", (GrassTensorAlgebraElement,), {"ring": self})

    def __hash__(self):
        return hash(("GrassTensorAlgebra", self._ID))

    def elem_sym(self, p, k):
        from schubmult import uncode
        set_of_keys = [(rc,) for rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1]*p), k)]
        return self.from_dict(dict.fromkeys(set_of_keys, S.One))

    def printing_term(self, key):
        return GrassTensorPrintingTerm(key)

    def _ensure_valid_rc_graph(self, rc: RCGraph, context: str = "") -> RCGraph:
        """Validate rc against RCGraph.all_rc_graphs(rc.perm, len(rc))."""
        valid_set = RCGraph.all_rc_graphs(rc.perm, len(rc))
        if rc not in valid_set:
            context_msg = f" ({context})" if context else ""
            raise ValueError(f"Invalid RCGraph produced{context_msg}: {rc}")
        return rc

    def key_to_rc_graph(self, key: tuple) -> RCGraph:
        """Evaluate a tensor key to an RCGraph using left-to-right squash_product."""
        if not isinstance(key, tuple):
            raise TypeError(f"Expected tuple key, got {type(key)}")
        if len(key) == 0:
            return RCGraph()

        # for rc in key:
        #     if not isinstance(rc, RCGraph):
        #         raise TypeError(f"Key factors must be RCGraph, got {type(rc)}")
        #     self._ensure_valid_rc_graph(rc, context="key_to_rc_graph input factor")

        acc = key[0]
        for factor in key[1:]:
            # Keep lengths compatible before squashing left to right.
            if len(acc) < len(factor):
                acc = acc.resize(len(factor))
            #     self._ensure_valid_rc_graph(acc, context="key_to_rc_graph resize accumulator")
            # elif len(factor) < len(acc):
            #     # factor = factor.resize(len(acc))
            #     # self._ensure_valid_rc_graph(factor, context="key_to_rc_graph resize factor")
            #     raise ValueError(f"Unexpected length decrease in key_to_rc_graph: acc length {len(acc)}, factor length {len(factor)}")

            acc = acc.squash_product(factor)
            # self._ensure_valid_rc_graph(acc, context="key_to_rc_graph squash step")

        return acc.normalize()

    def _normalize_key(self, key: tuple) -> tuple:
        """Normalize an RCGraph tuple to normal form."""
        _ensure_valid_key(key)
        result: list = [*key]
        # rows increasing
        # found = True
        # while found:
        #     found = False
        for i in range(1, len(result)):
            if len(result[i]) <= len(result[i - 1]):
                index = i
                while index > 0 and len(result[index]) < len(result[index - 1]):
                    result[index], result[index - 1] = result[index - 1], result[index]
                    index -= 1
                if index > 0 and len(result[index - 1]) == len(result[index]):
                    merged = result[index - 1].squash_product(result[index])
                    result[index - 1] = merged
                    del result[index]
                    index -= 1
                return self._normalize_key(tuple(result))
            if i < len(result) - 1:
                rc_left, rc_right = result[i].extend(1).squash_decomp()
                if rc_left.perm.inv !=0 and len(rc_left.perm.descents()) == 1 and rc_right.perm.inv != 0:
                    result[i] = rc_left.resize(_last_descent_size(rc_left))
                    result.insert(i + 1, rc_right)
                    return self._normalize_key(tuple(result))

        result = [rc for rc in result if rc.perm.inv != 0]  # drop identity factors
        # key should be increasing elem syms, and last can be arbitrary Grassmannian
        return tuple(result)

    def _mul_keys(self, left_key: tuple, right_key: tuple) -> tuple:
        new_key = []
        stck = list(reversed(left_key))
        stck2 = list(reversed(right_key))
        while len(stck) > 0 or len(stck2) > 0:
            if len(stck) == 0:
                new_key.append(stck2.pop())
            elif len(stck2) == 0:
                new_key.append(stck.pop())
            else:
                left_top = stck[-1]
                right_top = stck2[-1]
                if len(left_top) < len(right_top):
                    new_key.append(stck.pop())
                elif len(right_top) < len(left_top):
                    new_key.append(stck2.pop())
                else:
                    merged = left_top.squash_product(right_top)
                    stck.pop()
                    stck2.pop()
                    if merged.perm.inv != 0:
                        new_key.append(merged)

        return self._normalize_key(tuple(new_key))

    def from_dict(self, dct):
        accum = {}
        for key, coeff in dct.items():
            if coeff == S.Zero:
                continue
            nkey = self._normalize_key(key)
            accum[nkey] = accum.get(nkey, S.Zero) + coeff
        return self.dtype({k: v for k, v in accum.items() if v != S.Zero})

    def __call__(self, key):
        return self.from_dict({self._normalize_key(key): S.One})

    def from_rc_graph(self, rc: RCGraph):
        # """Create a single-key element from an arbitrary RC graph via normal-form factorization."""
        # if not isinstance(rc, RCGraph):
        #     raise TypeError(f"from_rc_graph expects RCGraph, got {type(rc)}")
        # key = self._factor_to_normal_tuple(rc)
        # assert rc == self.key_to_rc_graph(key), f"Factorization mismatch: {rc} != {self.key_to_rc_graph(key)}"
        # return self.from_dict({key: S.One})
        raise NotImplementedError("from_rc_graph is not implemented; use from_dict with explicit keys for now.")

    def mul(self, a, b):
        if not isinstance(b, GrassTensorAlgebraElement):
            return super().mul(a, b)
        result = {}
        for left_key, left_coeff in a.items():
            for right_key, right_coeff in b.items():
                key = self._mul_keys(left_key, right_key)
                result[key] = result.get(key, S.Zero) + left_coeff * right_coeff
        return self.from_dict(result)

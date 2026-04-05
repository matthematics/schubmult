from __future__ import annotations

from sympy import Tuple

from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.combinatorial.crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
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
    if not isinstance(key, (tuple, CrystalGraphTensor)):
        raise TypeError(f"Expected CrystalGraphTensor or tuple key, got {type(key)}")
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


class GrassTensorAlgebraElement(CrystalGraphRingElement):
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

    def dual_product(self, other):
        result = self.ring.zero
        for key1, coeff1 in self.items():
            for key2, coeff2 in other.items():
                result += coeff1 * coeff2 * self.ring.dual_product_on_basis(key1, key2)
        return result

class GrassTensorAlgebra(CrystalGraphRing):
    """Tensor-like algebra on tuples of full Grassmannian RC graphs.

    Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
    RC graph. Simplification rules
    """

    _id = 0
    zero_monom = CrystalGraphTensor()

    def __init__(self, domain=None):
        super().__init__(domain=domain)
        self._ID = GrassTensorAlgebra._id
        GrassTensorAlgebra._id += 1
        self.dtype = type("GrassTensorAlgebraElement", (GrassTensorAlgebraElement,), {"ring": self})

    def dual_product_on_basis(self, left_key, right_key):
        """Dual product to the coproduct_on_basis deconcatenation."""
        import itertools

        from schubmult.rings.combinatorial.rc_graph_ring import GrassRCGraphRing
        g = GrassRCGraphRing()
        result = self.one
        rev_left_key = reversed(left_key.factors)
        rev_right_key = reversed(right_key.factors)
        for left_factor, right_factor in itertools.zip_longest(rev_left_key, rev_right_key, fillvalue=RCGraph([])):
            result = result *self.from_dict({(k,): v for k, v in (g(left_factor) * g(right_factor)).items()})
        return self.from_dict(result)

    def coproduct_on_basis(self, key):
        if len(key) == 0:
            return (self@self)((key,key))
        result = (self@self).zero
        for i in range(len(key[-1]) + 1):
            left_key = CrystalGraphTensor(*[key[j].vertical_cut(i)[0] for j in range(len(key)) if len(key[j]) >= i])
            right_key = CrystalGraphTensor(*[key[j].vertical_cut(i)[1] for j in range(len(key)) if len(key[j]) >= i])
            result += self(left_key) @ self(right_key)
        return result

    def _ensure_cg(self, obj: CrystalGraphTensor | None) -> CrystalGraphTensor | None:
        """
        Ensure obj behaves like a CrystalGraph.
        """
        if obj is None:
            return None
        if isinstance(obj, CrystalGraphTensor):
            return obj
        if isinstance(obj, tuple):
            return CrystalGraphTensor(*obj)
        return obj

    def __hash__(self):
        return hash(("GrassTensorAlgebra", self._ID))

    def elem_sym(self, p, k):
        from schubmult import uncode
        set_of_keys = [CrystalGraphTensor(rc) for rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1]*p), k)]
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

    def key_to_rc_graph(self, key: CrystalGraphTensor | tuple) -> RCGraph:
        """Evaluate a tensor key to an RCGraph using left-to-right squash_product."""
        if not isinstance(key, (tuple, CrystalGraphTensor)):
            raise TypeError(f"Expected CrystalGraphTensor or tuple key, got {type(key)}")
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

    def _normalize_key(self, key: CrystalGraphTensor | tuple) -> CrystalGraphTensor:
        """Normalize an RCGraph tensor key to normal form."""
        if isinstance(key, tuple):
            key = CrystalGraphTensor(*key)
        _ensure_valid_key(key)
        result: list = [*key]
        # rows increasing
        # found = True
        # while found:
        #     found = False
        for i in range(1, len(result)):
            if len(result[i]) == len(result[i - 1]):
                merged = result[i - 1].squash_product(result[i])
                result[i - 1] = merged
                del result[i]
                return self._normalize_key(CrystalGraphTensor(*result))
            if len(result[i]) < len(result[i - 1]):
                result[i - 1], result[i] = result[i], result[i - 1]
                return self._normalize_key(CrystalGraphTensor(*result))
        result = [rc.normalize() for rc in result if rc.perm.inv != 0]  # drop identity factors
        # def _is_elem_sym(rc: RCGraph) -> bool:
        #     from schubmult import uncode
        #     return rc.perm == uncode([0] * (len(rc) - rc.perm.inv) + [1] * (rc.perm.inv))
        # if len(result) > 1:

            # for i in range(len(result) - 1):
            #     if not _is_elem_sym(result[i]):
            #         combined = result[i].resize(len(result[i + 1])).squash_product(result[i + 1])
            #         result_base, result_grass = combined.squash_decomp()
            #         result_base = result_base.normalize()
            #         if result_grass.perm.inv == 0:
            #             raise ValueError(f"Unexpected full squash to identity in _normalize_key: {result[i]} * {result[i+1]} = {combined}")
            #         if result_base.perm.inv == 0:
            #             result[i] = result_grass.normalize()
            #             del result[i + 1]
            #             return self._normalize_key(CrystalGraphTensor(*result))
            #         if result[i + 1] == result_grass:
            #             if _is_elem_sym(result_base):
            #                 result[i] = result_base
            #                 continue
            #         result[i] = result_base
            #         result[i + 1] = result_grass
            #         while not _is_elem_sym(result[i]):
            #             result_base, result_grass = result[i].squash_decomp()
            #             if result_base.perm.inv == 0:
            #                 try_length = len(result[i]) + 1
            #                 while try_length <= len(result[i+1]) and (result_base.perm.inv == 0 or result_grass.perm.inv == 0):
            #                     result_base, result_grass = result[i].resize(try_length).squash_decomp()
            #                     try_length += 1
            #                 #result_base, result_grass = result[i].resize(len(result[i]) + 1).squash_decomp()
            #             result[i] = result_grass
            #             if result_base.perm.inv == 0:
            #                 del result[i]
            #                 break
            #             result.insert(i, result_base.normalize())
            #         result = [rc.normalize() for rc in result if rc.perm.inv != 0]  # drop identity factors
            #         return self._normalize_key(CrystalGraphTensor(*result))


        # key should be increasing elem syms, and last can be arbitrary Grassmannian
        return CrystalGraphTensor(*result)

    def _mul_keys(self, left_key: tuple, right_key: tuple) -> tuple:
        # new_key = []
        # stck = list(reversed(left_key))
        # stck2 = list(reversed(right_key))
        # while len(stck) > 0 or len(stck2) > 0:
        #     if len(stck) == 0:
        #         new_key.append(stck2.pop())
        #     elif len(stck2) == 0:
        #         new_key.append(stck.pop())
        #     else:
        #         left_top = stck[-1]
        #         right_top = stck2[-1]
        #         if len(left_top) < len(right_top):
        #             new_key.append(stck.pop())
        #         elif len(right_top) < len(left_top):
        #             new_key.append(stck2.pop())
        #         else:
        #             merged = left_top.squash_product(right_top)
        #             stck.pop()
        #             stck2.pop()
        #             if merged.perm.inv != 0:
        #                 new_key.append(merged)
        new_key = left_key + right_key
        return self._normalize_key(new_key)

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

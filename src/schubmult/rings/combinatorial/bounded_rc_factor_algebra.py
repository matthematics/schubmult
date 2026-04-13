from __future__ import annotations

from functools import cache

from sympy import Tuple

from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.combinatorial.crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRing, RCGraphRingElement
from schubmult.rings.printing import PrintingTerm, TypedPrintingTerm
from schubmult.symbolic import S


def _is_full_grassmannian_rc(rc: RCGraph) -> bool:
    # try:
    #     attr = rc.is_full_grassmannian
    # except AttributeError:
    #     attr = None
    # if isinstance(attr, bool):
    #     return attr
    return rc.perm.inv == 0 or rc.perm.descents() == {len(rc) - 1}


def _descent_of_grass(rc: RCGraph) -> int:
    if rc.perm.inv == 0:
        return -1
    descs = rc.perm.descents()
    if len(descs) == 0:
        return -1
    return max(descs) + 1


def _last_descent_size(rc: RCGraph) -> int:
    """Return max descent + 1, or 0 for identity."""
    if rc.perm.inv == 0:
        return 0
    descs = rc.perm.descents()
    if len(descs) == 0:
        return 0
    return max(descs) + 1





class BoundedRCFactorPrintingTerm(PrintingTerm):
    is_commutative = False
    precedence = 50



    def __new__(cls, key):
        return BoundedRCFactorPrintingTerm.__xnew_cached__(cls, key)

    @staticmethod
    def __xnew__(_class, key):
        obj = PrintingTerm.__new__(_class, key, None, None)
        obj._key = key
        return obj

    @staticmethod
    def __xnew_cached__(_class, key):
        return BoundedRCFactorPrintingTerm.__xnew__(_class, key)

    def __hash__(self):
        return hash((self._key, "BoundedRCFactorPrintingTerm"))

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


class BoundedRCFactorAlgebraElement(CrystalGraphRingElement):
    """Element of BoundedRCFactorAlgebra: finite linear combinations of Grass tensors."""

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

    def divdiff_descs(self, *indexes):
        result = self.ring.zero
        for index in reversed(indexes):
            for key, coeff in self.items():
                result += coeff * self.ring.div_diff_desc_key(index, key)
        return result

    def divdiff_perm(self, perm):
        result = self.ring.zero
        for key, coeff in self.items():
            result += coeff * self.ring.div_diff_perm_key(perm, key)
        return result

class BoundedRCFactorAlgebra(CrystalGraphRing):
    """Tensor-like algebra on tuples of full Grassmannian RC graphs.

    Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
    RC graph. Simplification rules
    """

    _id = 0
    #zero_monom = BoundedRCFactorAlgebra._key((),0)

    @property
    def make_key(self):
        return BoundedRCFactorAlgebra._key

    class _key(CrystalGraphTensor):
        def __init__(self, the_tup, size):
            self._size = size
            super().__init__(*the_tup)

        @property
        def size(self):
            return self._size

        def __hash__(self):
            return hash((self._size, self.factors))

    def __init__(self, domain=None):
        super().__init__(domain=domain)
        self.zero_monom = self.make_key((),0)
        self._ID = BoundedRCFactorAlgebra._id
        BoundedRCFactorAlgebra._id += 1

        self.dtype = type("BoundedRCFactorAlgebraElement", (BoundedRCFactorAlgebraElement,), {"ring": self})

    @cache
    def _schub_elem_cached(self, perm, size):
        #dct = RCGraph.full_CEM(perm, size, partition=tuple((~(perm.mul_dominant())).trimcode))
        dct = RCGraph.full_CEM(perm, size)
        elem = sum([self.from_tensor_dict(cem_dict, size=size) for _, cem_dict in dct.items()])
        return elem

    def div_diff_desc_key(self, index, key):
        # r = RCGraphRing()
        key = self._normalize_key(key)
        if not any(_descent_of_grass(rc) == index for rc in key):
            return self.zero
        new_keys = []
        for pos, rc in enumerate(key):
            if _descent_of_grass(rc) < index:
                continue
            if _descent_of_grass(rc) == index:
                for new_rc in rc.divdiff_desc(index):
                    if not _is_full_grassmannian_rc(new_rc):
                        raise ValueError(f"divdiff_desc produced non full Grassmannian RC graph: {new_rc} from {rc} at index {index} {key=} {key.size=}")
                    new_keys.append(tuple(new_rc.normalize() if the_pos == pos else key[the_pos] for the_pos in range(len(key))))
            if _descent_of_grass(rc) > index:
                new_new_keys = []
                for new_key in new_keys:
                    new_new_keys.append(tuple(rc.crystal_reflection(index) if the_pos == pos else new_key[the_pos] for the_pos in range(len(key))))
                new_keys = new_new_keys
        return sum([self(self.make_key(new_key, key.size)) for new_key in new_keys])

    def div_diff_perm_key(self, perm, key):
        result = self(key)
        while True:
            descs = perm.descents()
            if len(descs) == 0:
                break
            index = max(descs) + 1
            new_result = self.zero
            for the_key, coeff in result.items():
                if coeff != 0:
                    new_result += coeff * self.div_diff_desc_key(index, the_key)
            result = new_result
            perm = perm.swap(index - 1, index)
        return result


    def schub_elem(self, perm, size):
        from schubmult.combinatorics.permutation import Permutation
        return self._schub_elem_cached(Permutation(perm), size)

    def dual_product_on_basis(self, left_key, right_key):
        """Dual product to the coproduct_on_basis deconcatenation."""
        import itertools

        r =RCGraphRing()
        result = self.one
        rev_left_key = reversed(left_key.factors)
        rev_right_key = reversed(right_key.factors)
        for left_factor, right_factor in itertools.zip_longest(rev_left_key, rev_right_key, fillvalue=RCGraph([])):
            result = self.from_tensor_dict({(k.normalize(),): v for k, v in (r(left_factor) * r(right_factor)).items() if _is_full_grassmannian_rc(k.normalize())}, size=left_key.size+right_key.size) * result
        return self.from_dict(result)

    def coproduct_on_basis(self, key):
        if key.size == 0:
            return (self@self)((key,key))
        result = (self@self).zero
        for i in range(key.size + 1):
            left_key = []
            right_key = []
            for j in range(len(key)):
                left_factor, right_factor = key[j].resize(key.size).vertical_cut(i)
                left_key.append(left_factor.normalize())
                right_key.append(right_factor.normalize())
            result += self(self.make_key(left_key, i)) @ self(self.make_key(right_key, key.size - i))
        return result

    # def _ensure_cg(self, obj: CrystalGraphTensor | None) -> CrystalGraphTensor | None:
    #     """
    #     Ensure obj behaves like a CrystalGraph.
    #     """
    #     if obj is None:
    #         return None
    #     if isinstance(obj, CrystalGraphTensor):
    #         return obj
    #     if isinstance(obj, tuple):
    #         return CrystalGraphTensor(*obj)
    #     return obj

    def __hash__(self):
        return hash(("BoundedRCFactorAlgebra", self._ID))

    def elem_sym(self, p, k):
        from schubmult import uncode
        set_of_keys = [self.make_key((rc,),k) for rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1]*p), k)]
        return self.from_dict(dict.fromkeys(set_of_keys, S.One))

    def printing_term(self, key):
        return BoundedRCFactorPrintingTerm(key)

    def key_to_rc_graph(self, key) -> RCGraph:
        """Evaluate a tensor key to an RCGraph using left-to-right squash_product."""
        if not isinstance(key, self.make_key):
            raise TypeError(f"Expected CrystalGraphTensor or tuple key, got {type(key)}")
        if len(key) == 0:
            return RCGraph([]).resize(key.size)

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

        return acc.resize(key.size)

    def _ensure_valid_key(self, key):
        if not isinstance(key, (tuple, self.make_key)):
            raise TypeError(f"Expected key type or tuple key, got {type(key)} for {key=}")
        if isinstance(key, tuple):
            key = self.make_key(key[0], key[1])
        for rc in key:
            if not isinstance(rc, RCGraph):
                raise TypeError(f"Key factors must be RCGraph, got {type(rc)}")
            if len(rc) > key.size:
                raise ValueError(f"Factor of size {len(rc)} in {key} which is bigger than {key.size}")
            if not _is_full_grassmannian_rc(rc):
                raise ValueError(f"Key factors must be full Grassmannian RC graphs, got {rc}")
        return key


    def _normalize_key(self, key):
        """Normalize an RCGraph tensor key to normal form."""
        # if any(len(k) > key.size for k in key):
        #     raise ValueError(f"Key factors too big in {key} with size {key.size}")

        key = self._ensure_valid_key(key)
        result: list = [*key]
        # rows increasing
        # found = True
        # while found:
        #     found = False
        for i in range(1, len(result)):
            if len(result[i]) == len(result[i - 1]):
                merged = result[i - 1].squash_product(result[i])
                result[i - 1] = merged
                #del result[i]
                result = result[:i] + result[i + 1:]
                return self._normalize_key(self.make_key(result, key.size))
            if len(result[i]) < len(result[i - 1]):
                result[i - 1], result[i] = result[i], result[i - 1]
                return self._normalize_key(self.make_key(result, key.size))
        result = [rc.normalize() for rc in result if rc.perm.inv != 0]  # drop identity factors

        #the_rc = self.key_to_rc_graph(self.make_key(result, key.size))
        for i in range(len(result) - 1, -1, -1):
            if not _is_full_grassmannian_rc(result[i]):
                raise ValueError(f"Non full Grassmannian factor in normalized key: {result[i]} in key {key}")
            if len(result[i].perm) - 1 > len(result[i]) and (len(result[i]) < key.size):
                max_len = min(len(result[i].perm), key.size)
                rc = result[i].resize(max_len)
                rc_base, rc_grass = rc.squash_decomp()
                rc_base = rc_base.normalize()
                rc_grass = rc_grass.normalize()
                if rc_base.perm.inv != 0 and rc_grass.perm.inv != 0 and _is_full_grassmannian_rc(rc_base):
                    result[i] = rc_base.normalize()
                    result.insert(i + 1, rc_grass.normalize())
                    return self._normalize_key(self.make_key(result, key.size))

        return self.make_key(result, key.size)

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
        new_key = self.make_key((*left_key, *right_key), max(left_key.size, right_key.size))
        return self._normalize_key(new_key)

    def from_dict(self, dct):
        accum = {}
        for key, coeff in dct.items():
            if coeff == S.Zero:
                continue
            nkey = self._normalize_key(key)
            accum[nkey] = accum.get(nkey, S.Zero) + coeff
        return self.dtype({k: v for k, v in accum.items() if v != S.Zero})

    def from_tensor_dict(self, dct, size):
        accum = {}
        for key, coeff in dct.items():
            if coeff == S.Zero:
                continue
            nkey = self._normalize_key(self.make_key(key, size))
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
        if not isinstance(b, BoundedRCFactorAlgebraElement):
            return super().mul(a, b)
        result = {}
        for left_key, left_coeff in a.items():
            for right_key, right_coeff in b.items():
                key = self._mul_keys(left_key, right_key)
                result[key] = result.get(key, S.Zero) + left_coeff * right_coeff
        return self.from_dict(result)

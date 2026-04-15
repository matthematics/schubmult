from __future__ import annotations

import logging
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

def _is_elem_sym_rc(rc: RCGraph) -> bool:
    return _is_full_grassmannian_rc(rc) and len(rc.perm) - 1 <= len(rc)


@cache
def _max_grass_elem_peel(rc):
    """Try to peel a wide grassmannian RC graph into a smaller piece + elementary symmetric piece."""
    from schubmult import uncode

    perm = rc.perm
    r = RCGraphRing()
    if _is_elem_sym_rc(rc):
        return None, None
    the_mully = r(rc) * (r.monomial(*([0] * (len(perm)))))
    length = len(rc) + 1
    max_found = -1
    the_min = None
    the_max = None
    for rcc in the_mully:
        new_rcc = RCGraph(rcc)
        rows = []
        maxxy = max(new_rcc.perm.descents()) + 1
        while True:
            if maxxy < 1:
                break
            if new_rcc.perm[maxxy - 1] > new_rcc.perm[maxxy]:
                new_rcc, row = new_rcc.exchange_property(maxxy, return_row=True)
                rows.append(row)
                # if len(rows) > max_found:
                #     break
                maxxy -= 1
            else:
                break
        if len(set(rows)) != len(rows):
            continue
        new_rcc = new_rcc.normalize()
        if len(rows) > max_found and _is_full_grassmannian_rc(new_rcc) and (len(new_rcc) < len(rc) or (len(rc) == len(new_rcc) and new_rcc.perm.bruhat_leq(rc.perm))):
            weight = [0] * length
            for row in rows:
                weight[row - 1] += 1
            elem_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (length - len(rows)) + [1] * len(rows)), length, weight=tuple(weight))))
            if new_rcc.resize(len(elem_rc)).squash_product(elem_rc).normalize() == rc:
                the_min, the_max = new_rcc, elem_rc
                max_found = elem_rc.inv
    return the_min, the_max


def _sort_and_merge(factors):
    """Sort factors by length (ascending, stable) and merge same-length adjacent via squash_product."""
    factors.sort(key=len)
    if len(factors) <= 1:
        return factors
    merged = [factors[0]]
    for rc in factors[1:]:
        if len(merged[-1]) == len(rc):
            m = merged[-1].squash_product(rc)
            if m.perm.inv == 0:
                merged.pop()
            else:
                merged[-1] = m
        else:
            merged.append(rc)

    merged = [m.normalize() for m in merged if m.perm.inv != 0]
    # while True:
    #     did_merge = False
    #     for i in range(len(merged) - 2, -1, -1):
    #         if i < len(merged) - 2:
    #             m = merged[i].resize(len(merged[i+1])).squash_product(merged[i + 1]).normalize()
    #             if _is_elem_sym_rc(m):
    #                 merged[i] = m
    #                 merged.pop(i + 1)
    #                 did_merge = True
    #                 break
    #         else:
    #             m = merged[i].resize(len(merged[i+1])).squash_product(merged[i + 1]).normalize()
    #             if _is_full_grassmannian_rc(m):
    #                 merged[i] = m
    #                 merged.pop(i + 1)
    #                 did_merge = True
    #                 break
    #     if not did_merge:
    #         break

    return merged


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
    # zero_monom = BoundedRCFactorAlgebra._key((),0)

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
        self.zero_monom = self.make_key((), 0)
        self._ID = BoundedRCFactorAlgebra._id
        BoundedRCFactorAlgebra._id += 1

        self.dtype = type("BoundedRCFactorAlgebraElement", (BoundedRCFactorAlgebraElement,), {"ring": self})

    @cache
    def _schub_elem_cached(self, perm, size, partition):
        if partition is None:
            partition = tuple((~(perm.mul_dominant())).trimcode)
        dct = RCGraph.full_CEM(perm, size, partition=partition)
        # dct = RCGraph.full_CEM(perm, size)
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

    def schub_elem(self, perm, size, partition = None):
        from schubmult.combinatorics.permutation import Permutation
        if partition is not None:
            partition = tuple(partition)
        return self._schub_elem_cached(Permutation(perm), size, partition=partition)

    def dual_product_on_basis(self, left_key, right_key):
        """Dual product to the coproduct_on_basis deconcatenation."""
        # import itertools

        # r =RCGraphRing()
        # result = self.one
        # # rev_left_key = reversed(left_key.factors)
        # # rev_right_key = reversed(right_key.factors)
        # index1 = 1
        # index2 = 1
        # the_index = 1
        # result = self.one
        # size = left_key.size + right_key.size
        # while the_index <= left_key.size + right_key.size:
        #     left_factor = left_key[index1 - 1] if (index1 <= len(left_key) and len(left_key[index1 - 1]) == the_index) else None
        #     right_factor = right_key[index2 - 1] if (index2 <= len(right_key) and len(right_key[index2 - 1]) == the_index) else None
        #     new_result = self.zero
        #     if left_factor is None and right_factor is None:
        #         the_index += 1
        #         new_result = result
        #     elif left_factor is None:
        #         new_result += result * self(((left_factor,),size))
        #         index1 += 1
        #         the_index += 1
        #     elif right_factor is None:
        #         new_result += result * self(((right_factor,),size))
        #         index2 += 1
        #         the_index += 1
        #     else:
        #         product = r(left_factor) * r(right_factor)
        #         new_result = self.zero
        #         for rc, coeff in product.items():
        #             if _is_full_grassmannian_rc(rc):
        #                 new_result += result * self(((rc,), size))
        #         index1 += 1
        #         index2 += 1
        #         the_index += 1
        #     result = new_result
        # return self.from_dict({(k, size): v for k, v in result.items()})

    def coproduct_on_basis(self, key):
        if key.size == 0:
            return (self @ self)((key, key))
        result = (self @ self).zero
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

        set_of_keys = [self.make_key((rc,), k) for rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k)]
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
                #logging.warning(f"Factor of size {len(rc)} in {key} which is bigger than {key.size}. Proceeding anyway")
                # rc = rc.resize(key.size)
            if not _is_full_grassmannian_rc(rc):
                raise ValueError(f"Key factors must be full Grassmannian RC graphs, got {rc}")
        return key

    def _normalize_key(self, key):
        """Normalize an RCGraph tensor key to normal form."""
        key = self._ensure_valid_key(key)
        size = key.size
        factors = list(key)

        while True:
            # Phase 1: sort by length and merge same-length adjacent factors
            factors = _sort_and_merge(factors)
            # Phase 2: normalize individual RCGraphs and strip identities
            factors = [rc.normalize() for rc in factors if rc.perm.inv != 0]

            # Phase 3: peel wide factors (perm wider than row count)
            peeled = False
            for i in range(len(factors) - 1):
                if not _is_full_grassmannian_rc(factors[i]):
                    raise ValueError(f"Non full Grassmannian factor in normalized key: {factors[i]} in key {key}")
                # if len(factors[i]) == size:
                #     break
                if len(factors[i].perm) - 1 > len(factors[i]) and (len(factors[i]) < size):
                    max_rc, elem_rc = _max_grass_elem_peel(factors[i])
                    if max_rc is not None and elem_rc is not None and elem_rc.perm.inv > 0 and max_rc.perm.inv > 0:
                        factors[i] = max_rc
                        factors.insert(i + 1, elem_rc)
                        peeled = True
                        break
            if not peeled:
                break
            #return self._normalize_key(self.make_key(factors, size))
        return self.make_key(factors, size)

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

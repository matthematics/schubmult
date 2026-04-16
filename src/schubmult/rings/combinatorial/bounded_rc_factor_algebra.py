from __future__ import annotations

import logging  # noqa: F401
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


@cache
def _max_grass_elem_peel(rc):
    """Try to peel a wide grassmannian RC graph into a smaller piece + elementary symmetric piece."""
    from schubmult import uncode

    # if len(rc) == 1:
    #     # special case
    #     if rc.perm.inv <= 1:
    #         return rc, RCGraph([])
    #     return RCGraph([rc[0][1:]]), RCGraph([(2,),()])
    perm = rc.perm
    r = RCGraphRing()
    # if rc.perm.inv == 1:
    if len(rc.perm) - 1 <= len(rc):
        return rc, RCGraph([])
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
                new_rcc_test, row = new_rcc.exchange_property(maxxy, return_row=True)
                if row in rows:
                    break
                new_rcc = new_rcc_test
                rows.append(row)
                maxxy -= 1
            else:
                break
        new_rcc = new_rcc.normalize()
        if len(rows) > max_found and _is_full_grassmannian_rc(new_rcc) and len(new_rcc.inv) < len(rc):
            #length = max(max(rows),min(len(rc.perm.trimcode) + 1, length))
            weight = [0] * length
            for row in rows:
                weight[row - 1] += 1
            elem_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (length - len(rows)) + [1] * len(rows)), length, weight=tuple(weight))))
            # if new_rcc.resize(len(elem_rc)).squash_product(elem_rc).normalize() == rc:
            the_min, the_max = new_rcc, elem_rc
            max_found = elem_rc.inv
    # if the_min is None:
    #     return _max_grass_elem_peel(rc.extend(1))
    return the_min, the_max


def _sort_and_merge(factors):
    """Sort factors by length (ascending, stable) and merge same-length adjacent via squash_product."""
    factors.sort(key=len)
    if len(factors) <= 1:
        return factors
    merged = [factors[0]]
    for rc in factors[1:]:
        if merged and len(merged[-1]) == len(rc):
            m = merged[-1].squash_product(rc)
            if m.perm.inv == 0:
                merged.pop()
            else:
                merged[-1] = m
        else:
            merged.append(rc)
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

        def __eq__(self, other):
            return isinstance(other, BoundedRCFactorAlgebra._key) and self.size == other.size and self.factors == other.factors

    def __init__(self, domain=None, post_normalize=False):
        super().__init__(domain=domain)
        self.post_normalize = post_normalize
        self.zero_monom = self.make_key((), 0)
        self._ID = BoundedRCFactorAlgebra._id
        BoundedRCFactorAlgebra._id += 1

        self.dtype = type("BoundedRCFactorAlgebraElement", (BoundedRCFactorAlgebraElement,), {"ring": self})

    @cache
    def _schub_elem_cached(self, perm, size, partition):
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

    def schub_elem(self, perm, size, partition=None):
        from schubmult.combinatorics.permutation import Permutation
        if size >= len(perm) - 1 and partition is None:
            partition = tuple(range(len(perm) - 1, 0, -1))
        elif partition is None:
            #partition = tuple((~(perm.strict_mul_dominant(size))).trimcode)
            partition = tuple((~(perm.strict_mul_dominant(size))).trimcode)

        return self._schub_elem_cached(Permutation(perm), size, partition)

    def dual_product_on_basis(self, left_key, right_key):
        """Dual product to the coproduct_on_basis deconcatenation."""
        r = RCGraphRing()
        total_rc = r(self.key_to_rc_graph(left_key)) * r(self.key_to_rc_graph(right_key))
        result = self.zero
        for big_rc, coeff in total_rc.items():
            partition = tuple((~(big_rc.perm.strict_mul_dominant(left_key.size + right_key.size))).trimcode)
            graphs = RCGraph.full_CEM(big_rc.perm, len(big_rc), partition=partition)
            cem = self.from_tensor_dict(graphs.get(big_rc, {}), size=len(big_rc))
            if cem.almosteq(self.zero):
                continue
            for cem_key, cem_coeff in cem.items():
                tensor_cut = [rc.vertical_cut(min(len(rc),left_key.size)) for rc in cem_key]
                top_cut_elem = self(self.make_key((rc1 for (rc1, rc2) in tensor_cut), size=left_key.size))
                bottom_cut_elem = self(self.make_key((rc2 for (rc1, rc2) in tensor_cut), size=right_key.size))
                top_cut = next(iter(top_cut_elem))
                bottom_cut = next(iter(bottom_cut_elem))
                if top_cut == left_key and bottom_cut == right_key:
                    result += self(cem_key)
        return result

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
                # logging.warning(f"Factor of size {len(rc)} in {key} which is bigger than {key.size}. Proceeding anyway")
                # rc = rc.resize(key.size)
            if not _is_full_grassmannian_rc(rc):
                raise ValueError(f"Key factors must be full Grassmannian RC graphs, got {rc}")
        return key

    def _normalize_key(self, key):
        """Normalize an RCGraph tensor key to normal form."""
        key = self._ensure_valid_key(key)
        size = key.size
        #hw_key, raise_seq = key.to_highest_weight()
        factors = list(key)
        #key_hw, raise_seq = key.to_highest_weight()
        #factors = list(key_hw)
        while True:
            # Phase 1: sort by length and merge same-length adjacent factors
            new_factors = _sort_and_merge(list(factors))
            #hw_key, raise_seq = self.make_key(tuple(new_factors), size).to_highest_weight()
            # Phase 2: normalize individual RCGraphs and strip identities
            factors = [rc.normalize() for rc in new_factors if rc.perm.inv != 0]

            peeled = False
            for i in range(len(factors)):
                if not _is_full_grassmannian_rc(factors[i]):
                    raise ValueError(f"Non full Grassmannian factor in normalized key: {factors[i]} in key {key}")
                factor = factors[i]
                if len(factor) < size:
                    max_rc, elem_rc = factor.resize(len(factor) + 1).squash_decomp()
                    max_rc = max_rc.normalize()
                    elem_rc = elem_rc.normalize()
                    if max_rc.perm.inv > 0 and elem_rc.perm.inv > 0 and len(elem_rc) <= size and _is_full_grassmannian_rc(max_rc) and _is_full_grassmannian_rc(elem_rc):
                        factors[i] = max_rc
                        factors.insert(i + 1, elem_rc)
                        peeled = True
                        break
        #         # if len(factors[i].perm) - 1 > len(factors[i]) and len(factors[i]) < size:
        #         #     factor = factors[i].resize(len(factors[i]) + 1)
        #         # else:
        #         if len(factors[i]) == size:
        #             continue
        #         if len(factors[i].perm) - 1 > len(factors[i]):
        #             factor = factors[i]
        #             max_rc, elem_rc = factor.squash_decomp()#_max_grass_elem_peel(factors[i])
        #             max_rc = max_rc.normalize()
        #             elem_rc = elem_rc.normalize()
        #             # if elem_rc.perm.inv > 0:
        #             #     raise ValueError(f"Peeling produced non-identity elem_rc: {elem_rc} from {factors[i]} in key {key}")
        #             if max_rc is not None and elem_rc is not None and elem_rc.perm.inv > 0 and max_rc.perm.inv > 0 and len(elem_rc) <= size and _is_full_grassmannian_rc(max_rc) and _is_full_grassmannian_rc(elem_rc):
        #                 factors[i] = max_rc
        #                 factors.insert(i + 1, elem_rc)
        #                 peeled = True
        #                 break
        #             # factor = factor.extend(1)
        #             # max_rc, elem_rc = _max_grass_elem_peel(factor)
        #             # max_rc = max_rc.normalize()
        #             # elem_rc = elem_rc.normalize()
        #             # # if elem_rc.perm.inv > 0:
        #             # #     raise ValueError(f"Peeling produced non-identity elem_rc: {elem_rc} from {factors[i]} in key {key}")
        #             # if max_rc is not None and elem_rc is not None and elem_rc.perm.inv > 0 and max_rc.perm.inv > 0 and len(elem_rc) <= size and _is_full_grassmannian_rc(max_rc) and _is_full_grassmannian_rc(elem_rc):
        #             #     factors[i] = max_rc
        #             #     factors.insert(i + 1, elem_rc)
        #             #     peeled = True
        #             #     break
        #         # if len(set(factors[i].perm.trimcode)) > 2:
        #         #     raise ValueError(f"Factor with more than 2 distinct row lengths in normalized key: {factors[i]} in key {key}")
        #     factors = self.make_key(self.make_key(tuple(factors), size).reverse_raise_seq(raise_seq), size)
            if not peeled:
                break

        # # for fact in factors:
        # #     if len(fact) < len(factors[-1]) and len(fact.perm) - 1 > len(fact):
        # #         raise ValueError(f"Factor {fact} in key {factors} is not wide enough to fit its permutation")
        #factors = _sort_and_merge(list(factors))
        return self.make_key(factors, size)
        #return self.make_key(tensor.reverse_raise_seq(raise_seq), size)

    def _mul_keys(self, left_key: tuple, right_key: tuple) -> tuple:
        if left_key.size != right_key.size:
            raise ValueError(f"Cannot multiply keys of different sizes: {left_key.size} vs {right_key.size}")
        new_key = self.make_key((*left_key, *right_key), left_key.size)
        return self._normalize_key(new_key)

    _post_normalizing = False

    def _post_normalize_dict(self, accum):
        """Expand keys whose factors are too wide (perm doesn't fit in row count) into sums via schub_elem.

        Uses a re-entrancy guard: if we're already inside a post-normalize call,
        just return the raw element to avoid infinite mutual recursion with
        mul -> _post_normalize_dict -> schub_elem -> from_tensor_dict -> _post_normalize_dict.
        """
        if self._post_normalizing or not self.post_normalize:
            return self.dtype({k: v for k, v in accum.items() if v != S.Zero})

        clean = {}
        to_expand = {}
        for key, coeff in accum.items():
            if coeff == S.Zero:
                continue
            bad = [i for i, factor in enumerate(key) if len(factor.perm) - 1 > len(factor) and len(factor) < key.size]
            if not bad:
                clean[key] = clean.get(key, S.Zero) + coeff
            else:
                to_expand[key] = coeff

        if not to_expand:
            return self.dtype({k: v for k, v in clean.items() if v != S.Zero})

        self._post_normalizing = False
        try:
            result = self.dtype({k: v for k, v in clean.items() if v != S.Zero})
            for key, coeff in to_expand.items():
                term = coeff * self.dtype({self.make_key((), key.size): S.One})
                for i in range(len(key)):
                    bad = len(key[i].perm) - 1 > len(key[i]) and len(key[i]) < key.size
                    if bad:
                        term = self.mul(term, self.schub_elem(key[i].perm, key.size))
                    else:
                        factor_elem = self.dtype({self.make_key((key[i],), key.size): S.One})
                        term = self.mul(term, factor_elem)
                for k, v in term.items():
                    result[k] = result.get(k, S.Zero) + v
        finally:
            self._post_normalizing = False

        # Re-check: the expansion may have produced new bad keys (one pass only)
        still_bad = any(
            len(factor.perm) - 1 > len(factor) and len(factor) < key.size
            for key in result
            for factor in key
            if result.get(key, S.Zero) != S.Zero
        )
        if still_bad:
            return self._post_normalize_dict(dict(result))
        return result

    def from_dict(self, dct):
        accum = {}
        for key, coeff in dct.items():
            if coeff == S.Zero:
                continue
            nkey = self._normalize_key(key)
            accum[nkey] = accum.get(nkey, S.Zero) + coeff
        return self._post_normalize_dict(accum)

    def from_tensor_dict(self, dct, size):
        accum = {}
        for key, coeff in dct.items():
            if coeff == S.Zero:
                continue
            nkey = self._normalize_key(self.make_key(key, size))
            accum[nkey] = accum.get(nkey, S.Zero) + coeff
        return self._post_normalize_dict(accum)

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
        accum = {}
        for left_key, left_coeff in a.items():
            for right_key, right_coeff in b.items():
                key = self._mul_keys(left_key, right_key)
                accum[key] = accum.get(key, S.Zero) + left_coeff * right_coeff
        return self._post_normalize_dict(accum)

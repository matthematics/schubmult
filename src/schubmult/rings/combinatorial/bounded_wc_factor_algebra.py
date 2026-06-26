from __future__ import annotations

import logging  # noqa: F401
from functools import cache

from sympy import Tuple

from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.combinatorial.crystal_graph_ring import CrystalGraphRing, CrystalGraphRingElement
from schubmult.rings.combinatorial.rc_graph_ring import RCGraphRing
from schubmult.rings.printing import PrintingTerm, TypedPrintingTerm
from schubmult.symbolic import S


def _tensor_to_rcs(weight_tensor, descents):
    from schubmult import uncode
    rcs = []
    for i, desc in enumerate(descents):
        wt = len(weight_tensor[i])
        weight = [0] * desc
        for w in weight_tensor[i]:
            weight[w - 1] = 1
        elem_rc = next(iter(RCGraph.all_wc_graphs(uncode([0]*(desc - wt) + [1] * wt), desc, weight=tuple(weight))))
        rcs.append(elem_rc)
    return CrystalGraphTensor(*rcs)

def _all_tensors(weights, descents):
    from schubmult.utils.schub_lib import hw_elementary_tensors
    hw_tensors = hw_elementary_tensors(weights, descents)
    for hw_tensor_weight in hw_tensors:
        hw_tensor = _tensor_to_rcs(hw_tensor_weight, descents)
        yield from hw_tensor.full_crystal

def _squash_it_up(tup):
    ret = RCGraph([()])
    for rc in tup:
        ret = ret.resize(len(rc)).squash_product(rc)
    return ret

def _elem_factor_from_rc(rc):
    from schubmult.combinatorics.permutation import Permutation
    if rc.perm.inv == 0:
        return {}
    n = len(rc.perm)
    weight = tuple(reversed([n - 1 - j - w for j, w in enumerate((rc.perm * Permutation.w0(n)).pad_code(n - 1))]))
    hw_rc, raise_seq = rc.to_highest_weight()
    good_tensor = None
    for tensor in _all_tensors(weight, tuple(range(1,n))):
        if _squash_it_up(tensor).resize(len(hw_rc)) == hw_rc:
            good_tensor = tensor
            break
    if good_tensor is None:
        raise ValueError(f"Could not find good tensor for WC graph {rc}")
    to_lower = CrystalGraphTensor(*good_tensor).reverse_raise_seq(raise_seq)
    return to_lower


def _is_full_grassmannian_rc(rc: RCGraph) -> bool:
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




def _build_elem_from_key(key):
    build_elem_key = []
    key_index = 0
    i = 1
    size = key.size
    while key_index < len(key):
        # a = key[i - 1]
        # if a != 0:
        #     build_elem_key.append(self.elem_sym(a, i, size=size))
        rc = key[key_index]
        if len(rc) == min(i, size):
            build_elem_key.append(rc.perm.inv)
            key_index += 1
        else:
            build_elem_key.append(0)
        i += 1
    if key_index != len(key):
        raise ValueError(f"Did not consume all of key {key} when building  {size} {key_index=} {len(key)=} {key=}")
    if len(build_elem_key) < size:
        build_elem_key.extend([0] * (size - len(build_elem_key)))
    elem_key = (tuple(build_elem_key[:size] + sorted(build_elem_key[size:], reverse=True)), size)
    return elem_key


def _build_schur_elem_from_key(key):
    build_elem_key = []
    key_index = 0
    i = 1
    size = key.size
    while key_index < len(key) and i < size:
        # a = key[i - 1]
        # if a != 0:
        #     build_elem_key.append(self.elem_sym(a, i, size=size))
        rc = key[key_index]
        if len(rc) == i:
            build_elem_key.append(rc.perm.inv)
            key_index += 1
        else:
            build_elem_key.append(0)
        i += 1
    # if key_index != len(key):
    #     raise ValueError(f"Did not consume all of key {key} when building  {size} {key_index=} {len(key)=} {key=}")
    if len(build_elem_key) < size - 1:
        build_elem_key.extend([0] * (size - 1 - len(build_elem_key)))
    partition = ()
    if key_index < len(key):
        partition = tuple(key[-1].perm.trimcode)
    else:
        partition = (0,) * size
    elem_key = (tuple(build_elem_key), partition)
    return elem_key




class BoundedWCFactorPrintingTerm(PrintingTerm):
    is_commutative = False
    precedence = 50

    def __new__(cls, key):
        return BoundedWCFactorPrintingTerm.__xnew_cached__(cls, key)

    @staticmethod
    def __xnew__(_class, key):
        obj = PrintingTerm.__new__(_class, key, None, None)
        obj._key = key
        return obj

    @staticmethod
    def __xnew_cached__(_class, key):
        return BoundedWCFactorPrintingTerm.__xnew__(_class, key)

    def __hash__(self):
        return hash((self._key, "BoundedWCFactorPrintingTerm"))

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


class BoundedWCFactorAlgebraElement(CrystalGraphRingElement):
    """Element of BoundedWCFactorAlgebra: finite linear combinations of Grass tensors."""

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else self[k] * self.ring.printing_term(k) for k in self.keys()]

    def to_rc_graph_ring_element(self):
        from schubmult import Gx
        from schubmult.combinatorics.rc_graph import RCGraph
        r = RCGraphRing()
        result = r.zero
        beta = Gx._beta
        for key, coeff in self.items():
            if coeff == 0:
                continue
            rc = self.ring.key_to_wc_graph(key)
            bagoinv = sum([len(wc.perm_word) - wc.perm.inv for wc in key])
            result += coeff * (beta**(bagoinv - len(rc.perm_word) + rc.perm.inv)) * r(RCGraph(rc))
        return result


class BoundedWCFactorAlgebra(CrystalGraphRing):
    """Tensor-like algebra on tuples of full Grassmannian WC graphs.

    Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
    WC graph. Simplification rules
    """

    _id = 0
    # zero_monom = BoundedWCFactorAlgebra._key((),0)

    @property
    def make_key(self):
        return BoundedWCFactorAlgebra._key

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
            return isinstance(other, BoundedWCFactorAlgebra._key) and self.size == other.size and self.factors == other.factors

    def __init__(self, domain=None, post_normalize=False):
        super().__init__(domain=domain)
        self.post_normalize = post_normalize
        self.zero_monom = self.make_key((), 0)
        self._ID = BoundedWCFactorAlgebra._id
        BoundedWCFactorAlgebra._id += 1

        self.dtype = type("BoundedWCFactorAlgebraElement", (BoundedWCFactorAlgebraElement,), {"ring": self})

    def from_rc_graph(self, rc, size):
        tensor = self._normalize_key(self.make_key(tuple(_elem_factor_from_rc(rc)), size))
        return self.from_dict({tensor: 1})

    def from_rc_graph_ring_element(self, elem, size):
        result = self.zero
        for rc, coeff in elem.items():
            result += coeff * self.from_rc_graph(rc, size)
        return result

    def from_CEM_rep(self, the_cem, size):
        from sympy import Add, Mul, Pow, expand, sympify

        # from schubmult import uncode

        terms = []
        the_cem = expand(sympify(the_cem))
        for the_term in Add.make_args(the_cem):
            coeff, rest = the_term.as_coeff_Mul()
            if rest == 1:
                terms.append((coeff, ()))
            else:
                factors = []
                for factor in Mul.make_args(rest):
                    exponent = 1
                    if isinstance(factor, Pow):
                        base, exponent = factor.as_base_exp()
                    else:
                        base = factor
                    factors.extend([base] * exponent)
                terms.append((coeff, sorted(factors, key=lambda x: (x.numvars, -x.degree))))

        the_result = self.zero
        for ti, term in enumerate(terms):
            base = self(self.make_key((), size))

            for fi, factor in enumerate(term[1]):
                p = factor.degree
                k = factor.numvars
                factor_result = self.elem_sym(p, k, size=size)

                base *= factor_result

            the_result += term[0] * base
        return the_result

    @cache
    def _schub_elem_cached(self, perm, size, partition):
        dct = RCGraph.full_CEM(perm, size, partition=partition)
        # dct = RCGraph.full_CEM(perm, size)
        elem = self.zero
        for _, cem_dict in dct.items():
            for key, coeff in cem_dict.items():
                new_key = self.make_key(key, size)
                # one_key = self.make_key((), size)
                # if not self._check_in_coprod(new_key, one_key, new_key):
                #     print(f"Key {new_key} from CEM of {perm} at size {size} failed coproduct check. Skipping.")
                #     continue
                elem += coeff * self(new_key)
        return elem


    def full_schub_elem(self, perm, size):
        # def cem_schub_schur_decomp(perm, n):
        #     from sympy import Add, Mul, Pow, expand, sympify

        #     from schubmult import Sx, uncode
        #     result = Sx.zero @ Sx.zero
        #     cd = (perm.strict_mul_dominant(n)).trimcode
        #     if any(a < n for a in cd):
        #         toadd = min(n - a for a in cd if a < n)
        #         cd = [a + toadd for a in cd]
        #     domperm = uncode(cd)
        #     reppy = sympify(expand(Sx(perm).cem_rep(mumu=~domperm, elem_func=Sx.symbol_elem_func), func=False))
        #     for arg in Add.make_args(reppy):
        #         coeff, schur_part = arg.as_coeff_Mul()
        #         part1 = Sx.one
        #         part2 = Sx.one
        #         for elem_arg in Mul.make_args(schur_part):
        #             if isinstance(elem_arg, Pow):
        #                 base, exp = elem_arg.as_base_exp()
        #             else:
        #                 base = elem_arg
        #                 exp = 1
        #             for _ in range(exp):
        #                 if base.numvars < n:
        #                     part1 *= base
        #                 else:
        #                     part2 *= base
        #             # if elem_arg.numvars < n:
        #             #     part1 *= elem_arg
        #             # else:
        #             #     part2 *= elem_arg
        #         result += coeff * part1 @ part2
        #     return result
        # decomp = cem_schub_schur_decomp(perm, size)
        from schubmult.rings.polynomial_algebra import ElemSymPolyBasis, Schub
        schub_elem = Schub(perm, size).change_basis(ElemSymPolyBasis)
        schub = self.zero
        for (comp, _), coeff in schub_elem.items():
            term = self(self.make_key((), size))
            #grass_part = self.make_key((), size)
            for index, part in enumerate(comp, start=1):
                if index < size:
                    term *= self.elem_sym(part, index, size=size)
                else:
                    term *= self.elem_sym(part, size, size=size)
            schub += coeff * term
        return schub

    def schub_elem(self, perm, size, partition=None):
        from schubmult.combinatorics.permutation import Permutation

        if partition is None:
            # partition = tuple((~(perm.strict_mul_dominant(size))).trimcode)
            partition = tuple((~(perm.strict_mul_dominant(size))).trimcode)

        return self._schub_elem_cached(Permutation(perm), size, partition)

    def schub_elem2(self, perm, size):
        raise NotImplementedError("schub_elem2 is looping over the wrong elem syms")
        from schubmult import ASx, ElementaryBasis, FreeAlgebra, SchubertBasis

        representation = ASx(perm, size).change_basis(ElementaryBasis)
        Elem = FreeAlgebra(ElementaryBasis)
        result = self.zero
        for key, coeff in representation.items():
            the_term = self(self.make_key((), size=size))
            key_tuple = key[0]
            the_elem = Elem(*key).change_basis(SchubertBasis)
            the_coeff = the_elem.get((perm, size), 0)
            if the_coeff == 0:
                continue
            for i, a in enumerate(key_tuple, start=1):
                if a != 0:
                    the_term *= self.elem_sym(a, min(i, size), size=size)
            result += the_coeff * coeff * the_term
        return result


    def __hash__(self):
        return hash(("BoundedWCFactorAlgebra", self._ID))

    def elem_sym(self, p, k, size):
        import math

        from schubmult import Gx, uncode
        res = self.zero
        beta = Gx._beta
        for i in range(p, k + 1):
            set_of_keys = [self.make_key((rc,), size) for rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k)]
            res += sum([beta**(i - p) * math.comb(i - 1, p - 1) * self.from_dict(dict.fromkeys(set_of_keys, S.One))])
        return res


    def printing_term(self, key):
        return BoundedWCFactorPrintingTerm(key)

    def key_to_wc_graph(self, key) -> RCGraph:
        """Evaluate a tensor key to an RCGraph using left-to-right squash_product."""
        if not isinstance(key, self.make_key):
            raise TypeError(f"Expected CrystalGraphTensor or tuple key, got {type(key)}\nfor {key=}")
        if len(key) == 0:
            return RCGraph([]).resize(key.size)


        acc = key[0]
        for factor in key[1:]:
            # Keep lengths compatible before squashing left to right.
            if len(acc) < len(factor):
                acc = acc.resize(len(factor))

            acc = acc.squash_product(factor)

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
                raise ValueError(f"Key factors must be full Grassmannian WC graphs, got {rc} for {key=}")
        return key

    def _merge_elem_sym(self, normalized_key, elem_sym_rc):
        if not _is_full_grassmannian_rc(elem_sym_rc):
            raise ValueError(f"elem_sym_rc must be full Grassmannian, got {elem_sym_rc}")
        if elem_sym_rc.perm.inv == 0:
            return normalized_key
        size = normalized_key.size
        if len(elem_sym_rc) > size:
            raise ValueError(f"elem_sym_rc of size {len(elem_sym_rc)} cannot fit into key of size {size}")
        if len(normalized_key) == 0:
            return (elem_sym_rc,)
        if len(elem_sym_rc) == size:
            if len(normalized_key[-1]) == size:
                return (*normalized_key[:-1], normalized_key[-1].squash_product(elem_sym_rc))
            return (*normalized_key, elem_sym_rc)
        if not any(len(rc) == len(elem_sym_rc) for rc in normalized_key):
            index = min([i for i, rc in enumerate(normalized_key) if len(rc) > len(elem_sym_rc)], default=len(normalized_key))
            return [*normalized_key[:index], elem_sym_rc, *normalized_key[index:]]
        index = next(i for i, rc in enumerate(normalized_key) if len(rc) == len(elem_sym_rc))
        tensor = CrystalGraphTensor(*normalized_key[:index + 1], elem_sym_rc)
        hw, raise_seq = tensor.to_highest_weight()
        left_part, the_elem_sym = hw.factors[:-1], hw.factors[-1]
        squashed_elem = left_part[-1].squash_product(the_elem_sym)
        left_part = left_part[:-1]
        base, overflow = squashed_elem.resize(len(elem_sym_rc) + 1).squash_decomp()
        base = base.normalize()
        overflow = overflow.normalize()
        build_back_up = self.make_key(left_part, size)#, base, *normalized_key[index + 1 :]
        if not _is_full_grassmannian_rc(base) or len(base) - 1 > len(base):
            base0 = base
            addup_list = []
            while base0.perm.inv > 0 and (not _is_full_grassmannian_rc(base0) or len(base0.perm) - 1 > len(base0)):
                base0, overflow0 = base0.squash_decomp()
                base0 = base0.normalize()
                overflow0 = overflow0.normalize()
                addup_list = [overflow0, *addup_list]
            if base0.perm.inv > 0:
                addup_list = [base0, *addup_list]
            for rc in addup_list:
                build_back_up = self.make_key(self._merge_elem_sym(build_back_up, rc), size)
        else:
            build_back_up = self.make_key(self._merge_elem_sym(build_back_up, base), size)
        build_back_up = self.make_key(self._merge_elem_sym(build_back_up, overflow), size)
        back_to_key = self.make_key(build_back_up.reverse_raise_seq(raise_seq), size)
        build_back_up = back_to_key
        for rc in normalized_key[index + 1:]:
            build_back_up = self.make_key(self._merge_elem_sym(build_back_up, rc), size)
        ret_key = list(build_back_up)
        return ret_key




    def _sort_and_merge(self, key, legacy=True):
        """Sort factors by length (ascending, stable) and merge same-length adjacent via squash_product."""
        if not legacy:
            new_key = self.make_key((), key.size)
            #last_seen = -1
            for index, rc in enumerate(key):
                if len(rc) == 0 or rc.perm.inv == 0:
                    continue
                new_key = self.make_key(self._merge_elem_sym(new_key, rc), key.size)
                #last_seen = 0 if len(new_key) == 0 else len(new_key[-1])
            if tuple(sorted(new_key, key=len)) != tuple(new_key) or not all(len(rc.perm) - 1 == len(rc) for rc in new_key if len(rc) < key.size):
                raise ValueError(f"Sorting and merging failed to produce length-sorted key: {new_key} from {key}")
            return new_key
        # legacy
        factors = list(key)
        #factors.sort(key=lambda )
        if len(factors) <= 1:
            return factors
        merged = [factors[0].normalize()]
        for rc in factors[1:]:
            rc = rc.normalize()
            index = len(merged) - 1
            while index >= 0 and len(merged[index]) > len(rc):
                index -= 1
            if index >= 0:
                # if len(merged[index]) == len(rc):
                #     m = merged[index].squash_product(rc)
                #     merged[index] = m
                # else:
                merged.insert(index + 1, rc)
            else:
                merged.insert(0, rc)
        return merged

    def _check_normal_key(self, key):
        if tuple(sorted(key, key=len)) != tuple(key):
            raise ValueError(f"Key is not sorted by length: {key}")
        if not all(len(rc.perm) - 1 == len(rc) for rc in key if len(rc) < key.size):
            raise ValueError(f"Key has non-identity factor that is not full Grassmannian: {key}")
        if sorted({len(rc) for rc in key}) != [len(rc) for rc in key]:
            raise ValueError(f"Key has factors of different sizes: {key}")
        return key


    def _normalize_key(self, key, legacy=True):
        """Normalize an RCGraph tensor key to normal form."""
        key = self._ensure_valid_key(key)
        if len(key) == 0:
            return key
        size = key.size
        factors = list(key)
        while True:
            factors = [rc.normalize() for rc in factors if rc.perm.inv != 0]
            factors = self._sort_and_merge(list(factors), legacy)
            factors = [rc.normalize() for rc in factors if rc.perm.inv != 0]
            break

        return self.make_key(factors, size)

    def _mul_keys(self, left_key: tuple, right_key: tuple) -> tuple:
        size = left_key.size
        if left_key.size != right_key.size:
            size = max(left_key.size, right_key.size)
        new_key = self._normalize_key(self.make_key((*left_key, *right_key), size))
        return new_key

    _post_normalizing = False

    def _post_normalize_dict(self, accum):
        """Expand keys whose factors are too wide (perm doesn't fit in row count) into sums via schub_elem.

        Uses a re-entrancy guard: if we're already inside a post-normalize call,
        just return the raw element to avoid infinite mutual recursion with
        mul -> _post_normalize_dict -> schub_elem -> from_tensor_dict -> _post_normalize_dict.
        """
        if True:
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
        still_bad = any(len(factor.perm) - 1 > len(factor) and len(factor) < key.size for key in result for factor in key if result.get(key, S.Zero) != S.Zero)
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


    def mul(self, a, b):
        if isinstance(a, BoundedWCFactorAlgebraElement):
            if isinstance(b, BoundedWCFactorAlgebraElement):
                accum = self.zero
                for left_key, left_coeff in a.items():
                    for right_key, right_coeff in b.items():
                        key = self._mul_keys(left_key, right_key)
                        accum += left_coeff * right_coeff * self(key)
                return accum
            return self.from_dict({k: v * b for k, v in a.items()})
        return self.from_dict({k: v * a for k, v in b.items()})


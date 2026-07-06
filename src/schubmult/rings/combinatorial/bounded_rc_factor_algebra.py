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


def _tensor_to_rcs(weight_tensor, descents):
    from schubmult import uncode

    rcs = []
    for i, desc in enumerate(descents):
        wt = len(weight_tensor[i])
        weight = [0] * desc
        for w in weight_tensor[i]:
            weight[w - 1] = 1
        elem_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (desc - wt) + [1] * wt), desc, weight=tuple(weight))))
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


def _is_row_root(root, k):
    if root[0] < k and root[1] >= k:
        return True
    return False


@cache
def _elem_factor_from_rc(rc):
    from schubmult.combinatorics.permutation import Permutation
    from schubmult.utils.perm_utils import has_bruhat_descent

    if rc.perm.inv == 0:
        return ()
    if rc.is_elem_sym:
        return CrystalGraphTensor(rc.normalize())
    n = len(rc.perm)
    k = n - 1
    weight = tuple(reversed([n - 1 - j - w for j, w in enumerate((rc.perm * Permutation.w0(n)).pad_code(n - 1))]))
    # hw_rc, raise_seq = rc.to_highest_weight()
    # good_tensor = None
    # for tensor in _all_tensors(weight, tuple(range(1,n))):
    #     if _squash_it_up(tensor).resize(len(hw_rc)) == hw_rc:
    #         good_tensor = tensor
    #         break
    degree = weight[-1]
    if degree == 0:
        raise ValueError(f"This shouldn't happen {rc=} {degree=}")
    stack = [[[], rc]]
    while stack:
        spot_list, test_rc = stack.pop()
        if len(spot_list) == degree:
            weight = [0] * k
            for r, c in spot_list:
                weight[r - 1] = 1
            elem_sym_rc = next(iter(RCGraph.elem_sym_rcs(degree, k, weight=tuple(weight))))
            if test_rc.resize(k).squash_product(elem_sym_rc).resize(len(rc)) == rc:
                return CrystalGraphTensor(*_elem_factor_from_rc(test_rc), elem_sym_rc)
            continue
        for i in range(k):
            for j in range(k, len(test_rc.perm)):
                if has_bruhat_descent(test_rc.perm, i, j):
                    row, col = test_rc.loc_of_inversion(i + 1, j + 1)
                    if row in spot_list:
                        continue
                    new_rc = test_rc.toggle_ref_at(row, col)
                    stack.append([[*spot_list, (row, col)], new_rc])

    # for the_spots in itertools.combinations(range(1, len(hw_rc) + 1), degree):
    #     test_rc = hw_rc
    #     rows = the_spots
    #     good = True
    #     for spot in the_spots:
    #         column =

    #     weight = [0] * k
    #     for r in rows:
    #         weight[r - 1] = 1
    #     elem_sym_rc = next(iter(RCGraph.elem_sym_rcs(degree, k, weight=tuple(weight)))).resize(len(rc))
    #     if test_rc.resize(k).squash_product(elem_sym_rc).resize(len(hw_rc)) == hw_rc:
    #         tensor = CrystalGraphTensor(test_rc, elem_sym_rc).reverse_raise_seq(raise_seq)
    #         return CrystalGraphTensor(*_elem_factor_from_rc(tensor.factors[0]), tensor.factors[1])
    # if good_tensor is None:
    raise ValueError(f"Could not find good tensor for RC graph {rc=} {degree=} {n - 1=}")
    # to_lower = CrystalGraphTensor(*good_tensor).reverse_raise_seq(raise_seq)
    # return to_lower


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

    def _to_top_forest_weight_graph(self, rc):
        return next(
            iter(
                [
                    rc2
                    for rc2 in RCGraph.all_forest_rcs(rc.forest_weight, weight=rc.length_vector)
                    if rc2.forest_weight == rc.forest_weight and rc2.omega_invariant[1] == rc.omega_invariant[1] and rc2.perm.pad_code(len(rc2)) == rc2.forest_weight
                ],
            ),
        )

    def prune(self):
        """Merge terms whose evaluated RC graphs share forest/omega invariants.

        Terms are grouped by
        ``(rc.forest_weight, rc.omega_invariant[1])`` where
        ``rc = self.ring.key_to_rc_graph(key)``.
        One key per group is kept (first encountered), and coefficients are summed.
        """
        representative_for_signature = {}
        coeff_by_representative = {}

        for key, coeff in self.items():
            if coeff == 0:
                continue
            rc = self._to_top_forest_weight_graph(self.ring.key_to_rc_graph(key))
            signature = rc  # (rc.forest_weight, rc.length_vector, rc.omega_invariant[1])
            representative = representative_for_signature.get(signature)
            if representative is None:
                representative_for_signature[signature] = key
                coeff_by_representative[key] = coeff
            else:
                coeff_by_representative[representative] += coeff

        result = self.ring.zero
        for key, coeff in coeff_by_representative.items():
            if coeff != 0:
                result += coeff * self.ring(key)
        return result

    def dual_product(self, other):
        result = self.ring.zero
        for key1, coeff1 in self.items():
            for key2, coeff2 in other.items():
                result += coeff1 * coeff2 * self.ring.dual_product_on_basis(key1, key2)
        return result

    def divdiff_descs(self, *indexes):
        result = self
        for index in reversed(indexes):
            old_result = result
            result = self.ring.zero
            for key, coeff in old_result.items():
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

        elem = self.zero
        for _, cem_dict in dct.items():
            for key, coeff in cem_dict.items():
                new_key = self.make_key(key, size)
                elem += coeff * self(new_key)
        return elem

    def full_schub_elem(self, perm, size):
        from schubmult.rings.polynomial_algebra import ElemSymPolyBasis, Schub

        schub_elem = Schub(perm, size).change_basis(ElemSymPolyBasis)
        schub = self.zero
        for (comp, _), coeff in schub_elem.items():
            term = self(self.make_key((), size))
            # grass_part = self.make_key((), size)
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
                tensor_cut = [rc.vertical_cut(min(len(rc), left_key.size)) for rc in cem_key]
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

    def elem_sym(self, p, k, size):
        from schubmult import uncode

        # size = k if size is None else size
        set_of_keys = [self.make_key((rc,), size) for rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), k)]
        # result = self.zero
        # if coeffvars is not None:
        #     for key in set_of_keys:
        #         result += self(key)
        #         coords = [key[0].left_to_right_inversion_coords(i) for i in range(key[0].perm.inv)]
        #         for r in range(1, p + 1):
        #             for spots in itertools.combinations(coords, r):
        #                 the_rc = key[0]
        #                 subtract = 0
        #                 for row, col in reversed(coords):
        #                     if (row, col) in spots:
        #                         subtract += 1
        #                         the_rc = the_rc.toggle_ref_at(row, col)
        #                     else:
        #                         if subtract > 0:
        #                             the_rc = the_rc.toggle_ref_at(row, col)
        #                             the_rc = the_rc.toggle_ref_at(row, col + subtract)
        #                 prd = prod([-coeffvars[col - 1] for row, col in spots])
        #                 if prd != 0:
        #                     result += prd * self(self.make_key((the_rc.normalize(),), size))
        # return result
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

    # def key_to_rc_graph_with_coeff(self, key, coeff_genset) -> RCGraph:
    #     """Evaluate a tensor key to an RCGraph using left-to-right squash_product."""
    #     import itertools
    #     from sympy import prod
    #     if not isinstance(key, self.make_key):
    #         raise TypeError(f"Expected CrystalGraphTensor or tuple key, got {type(key)}")
    #     if len(key) == 0:
    #         return RCGraph([]).resize(key.size)

    #     # for rc in key:
    #     #     if not isinstance(rc, RCGraph):
    #     #         raise TypeError(f"Key factors must be RCGraph, got {type(rc)}")
    #     #     self._ensure_valid_rc_graph(rc, context="key_to_rc_graph input factor")
    #     acc = self(self.make_key((), key.size))
    #     for i in range(len(key)):
    #         coords = [key[i].left_to_right_inversion_coords(j) for j in range(key[i].perm.inv)]
    #         result = self.zero
    #         for r in range(1, key[i].perm.inv + 1):
    #             for spots in itertools.combinations(coords, r):
    #                 the_rc = key[i]
    #                 subtract = 0
    #                 for row, col in reversed(coords):
    #                     if (row, col) in spots:
    #                         subtract += 1
    #                         the_rc = the_rc.toggle_ref_at(row, col)
    #                     else:
    #                         if subtract > 0:
    #                             the_rc = the_rc.toggle_ref_at(row, col)
    #                             the_rc = the_rc.toggle_ref_at(row, col + subtract)
    #                 coeffvars = coeff_genset[key.size - len(key[i]):]
    #                 prd = prod([-coeffvars[col - 1] for row, col in spots])
    #                 if prd != 0:
    #                     result += prd * self(self.make_key((the_rc.normalize(),), key.size))
    #         acc *= result

    #     return acc.to_rc_graph_ring_element()

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
        tensor = CrystalGraphTensor(*normalized_key[: index + 1], elem_sym_rc)
        hw, raise_seq = tensor.to_highest_weight()
        left_part, the_elem_sym = hw.factors[:-1], hw.factors[-1]
        squashed_elem = left_part[-1].squash_product(the_elem_sym)
        left_part = left_part[:-1]
        base, overflow = squashed_elem.resize(len(elem_sym_rc) + 1).squash_decomp()
        base = base.normalize()
        overflow = overflow.normalize()
        build_back_up = self.make_key(left_part, size)  # , base, *normalized_key[index + 1 :]
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
        for rc in normalized_key[index + 1 :]:
            build_back_up = self.make_key(self._merge_elem_sym(build_back_up, rc), size)
        #     full_overflow = addup_list
        # else:
        #     full_overflow = [base, overflow]
        # raise ValueError(f"Unexpected non full Grassmannian base from merging elem_sym: {base} from {elem_sym_rc} and {normalized_key[index]} in key of size {size}")
        # base0 = base
        # #full_overflow = [overflow]

        # while len(overflow) <= size and base0.perm.inv > 0 and (not _is_full_grassmannian_rc(base0) or len(base0.perm) - 1 > len(base0)):
        #     base0, overflow = base0.squash_decomp()
        #     base0 = base0.normalize()
        #     overflow = overflow.normalize()
        #     # merge partial key
        #     full_overflow = [overflow, *full_overflow]
        #     if overflow.perm.inv == 0:
        #         raise ValueError(f"Unexpected identity overflow during merging elem_sym: {overflow} from {base} and {normalized_key[index]} in key of size {size}")
        # base = base0
        # if base.perm.inv == 0:
        #     raise ValueError(f"Unexpected identity base from merging elem_sym: {base} from {elem_sym_rc} and {normalized_key[index]} in key of size {size}")
        # if len(base.perm) - 1 > len(base):
        #     raise ValueError(f"Unexpected non-identity base from merging elem_sym: {base} from {elem_sym_rc} and {normalized_key[index]} in key of size {size}")
        # mix_index = min([i for i, rc in enumerate(normalized_key) if len(rc) > len(base)], default=len(normalized_key))
        # full_overflow = [base, *full_overflow]
        # ret_key = [*normalized_key[:index]]
        # #the_rest = []
        # for mixer in [base, *normalized_key[index+1:]]:
        #     overflow_index = max([i for i, rc in enumerate(full_overflow) if len(rc) <= len(mixer)], default=-1) + 1
        #     #if overflow_index != -1:
        #     full_overflow = [*full_overflow[:overflow_index], mixer, *full_overflow[overflow_index:]]
        #     # else:
        #     #     full_overflow = [mixer, *full_overflow]
        # #flag = True
        # while len(full_overflow) > 0:
        #     overflow = full_overflow.pop(0)
        #     ret_key = self._merge_elem_sym(self.make_key(ret_key, size), overflow)
        # if index < len(normalized_key) - 1:
        #     if len(normalized_key[index + 1]) == len(overflow):
        #         ret_key = self.make_key([*self._merge_elem_sym(ret_key, overflow)], size)
        #         #, size), normalized_key[index + 1]), *normalized_key[index + 2:]]
        #     return (*normalized_key[:index], base, overflow, *normalized_key[index + 1 :])
        # return (*normalized_key[:index], base, overflow)
        ret_key = list(build_back_up)
        return ret_key
        # (*ret_key, ) if index < len(normalized_key) - 1 else ret_key
        # elif index < len(normalized_key) - 1 and len(normalized_key[index + 1]) == len(overflow):
        #     return [*normalized_key[:index], base, overflow, *normalized_key[index + 1:]]

    def _sort_and_merge(self, key):
        """Sort factors by length (ascending, stable) and merge same-length adjacent via squash_product."""
        factors = list(key)
        if len(factors) <= 1:
            return factors
        merged = [factors[0].normalize()]
        for rc in factors[1:]:
            rc = rc.normalize()
            index = len(merged) - 1
            while index >= 0 and len(merged[index]) > len(rc):
                index -= 1
            if index >= 0:
                if len(merged[index]) == len(rc):
                    m = merged[index].squash_product(rc)
                    merged[index] = m
                else:
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

    def _normalize_key(self, key):
        """Normalize an RCGraph tensor key to normal form."""
        key = self._ensure_valid_key(key)
        if len(key) == 0:
            return key

        size = key.size
        factors = list(key)

        factors = [rc.normalize() for rc in factors if rc.perm.inv != 0]
        factors = self._sort_and_merge(list(factors))
        factors = [rc.normalize() for rc in factors if rc.perm.inv != 0]

        return self.make_key(factors, size)

    def _check_in_coprod(self, left_key, right_key, new_key):
        from schubmult import ASx, ElementaryBasis, FreeAlgebra, SchubertBasis, SchurElementaryBasis  # noqa: F401

        SchurElem = FreeAlgebra(SchurElementaryBasis)  # noqa: F841

        rc_result = self.key_to_rc_graph(new_key)
        if ASx(rc_result.perm, len(rc_result)).change_basis(SchurElementaryBasis).coproduct().get((_build_schur_elem_from_key(left_key), _build_schur_elem_from_key(right_key)), 0) == 0:
            return False
        return True

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
        return self.dtype({k: v for k, v in accum.items() if v != S.Zero})

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
        if isinstance(a, BoundedRCFactorAlgebraElement):
            if isinstance(b, BoundedRCFactorAlgebraElement):
                accum = self.zero
                for left_key, left_coeff in a.items():
                    for right_key, right_coeff in b.items():
                        key = self._mul_keys(left_key, right_key)
                        # if key is not None:
                        accum += left_coeff * right_coeff * self(key)
                return accum
            return self.from_dict({k: v * b for k, v in a.items()})
        return self.from_dict({k: v * a for k, v in b.items()})

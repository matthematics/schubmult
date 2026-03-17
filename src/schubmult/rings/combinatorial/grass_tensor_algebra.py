from __future__ import annotations

from sympy import Tuple

from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.rings.base_ring import BaseRing, BaseRingElement
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

    def to_rc_graph(self, require_single_term: bool = True) -> RCGraph:
        """Evaluate to a single RC graph by left-to-right squash_product.

        If `require_single_term` is True, this requires exactly one basis key with
        coefficient 1 and returns its RCGraph evaluation.
        """
        if require_single_term:
            if len(self) != 1:
                raise ValueError(f"Expected exactly one term, got {len(self)} terms.")
            key, coeff = next(iter(self.items()))
            if coeff != S.One:
                raise ValueError(f"Expected coefficient 1 for single-term evaluation, got {coeff}.")
            return self.ring.key_to_rc_graph(key)

        raise NotImplementedError("Linear combination -> single RCGraph is not canonical; use require_single_term=True.")


class GrassTensorAlgebra(BaseRing):
    """Tensor-like algebra on tuples of full Grassmannian RC graphs.

    Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
    RC graph in strict increasing-length normal form.
    """

    _id = 0
    zero_monom = ()

    def __init__(self, domain=None, product_mode: str = "smash", compare_modes: bool = False):
        super().__init__(domain=domain)
        self._ID = GrassTensorAlgebra._id
        GrassTensorAlgebra._id += 1
        self.dtype = type("GrassTensorAlgebraElement", (GrassTensorAlgebraElement,), {"ring": self})
        if product_mode not in {"legacy", "smash", "forced"}:
            raise ValueError(f"Unknown product_mode={product_mode}; expected 'legacy', 'smash', or 'forced'.")
        self.product_mode = product_mode
        self.compare_modes = compare_modes

    def __hash__(self):
        return hash(("GrassTensorAlgebra", self._ID))

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

        for rc in key:
            if not isinstance(rc, RCGraph):
                raise TypeError(f"Key factors must be RCGraph, got {type(rc)}")
            self._ensure_valid_rc_graph(rc, context="key_to_rc_graph input factor")

        acc = key[0]
        for factor in key[1:]:
            # Keep lengths compatible before squashing left to right.
            if len(acc) < len(factor):
                acc = acc.resize(len(factor))
                self._ensure_valid_rc_graph(acc, context="key_to_rc_graph resize accumulator")
            elif len(factor) < len(acc):
                # factor = factor.resize(len(acc))
                # self._ensure_valid_rc_graph(factor, context="key_to_rc_graph resize factor")
                raise ValueError(f"Unexpected length decrease in key_to_rc_graph: acc length {len(acc)}, factor length {len(factor)}")

            acc = acc.squash_product(factor)
            self._ensure_valid_rc_graph(acc, context="key_to_rc_graph squash step")

        return acc.normalize()

    def same_size_product(self, left_grass: RCGraph, right_grass: RCGraph) -> RCGraph:
        """Product of full Grassmannian RC graphs of the same size."""
        if len(left_grass) != len(right_grass):
            raise ValueError(f"same_size_product requires equal lengths, got {len(left_grass)} and {len(right_grass)}")
        self._ensure_valid_rc_graph(left_grass, context="same_size_product left")
        self._ensure_valid_rc_graph(right_grass, context="same_size_product right")
        result = left_grass.left_squash(right_grass)
        self._ensure_valid_rc_graph(result, context="same_size_product result")
        return result

    def push_factor_through(self, left_grass: RCGraph, key: tuple) -> tuple:
        """Push a full Grassmannian factor through all adjacent smaller factors.

        This evaluates the maximal initial prefix of `key` with size strictly
        less than `left_grass`, right-squashes that prefix, left-squashes with
        `left_grass`, then right-squashes with the full suffix, and finally
        canonically expands back to Grass factors.
        """
        self._ensure_valid_rc_graph(left_grass, context="push_factor_through left")
        233
        if left_grass.perm.inv == 0:
            return key

        assert left_grass.is_full_grass, "Left factor must be full Grassmannian for push_factor_through"

        if len(key) == 0:
            return (left_grass,)

        # if len(key) == 1:
        #     right_factor = key[0]
        #     assert right_factor.is_full_grass, "Right factor must be full Grassmannian for push_factor_through single factor case"
        #     if len(right_factor) < len(left_grass):
        #         right_factor = right_factor.resize(len(left_grass))
        #         self._ensure_valid_rc_graph(right_factor, context="push_factor_through single factor resize")
        #         result = left_grass.left_squash(right_factor)
        #     else:
        #         left_factor = left_grass.resize(len(right_factor))
        #         self._ensure_valid_rc_graph(left_factor, context="push_factor_through single factor resize")
        #         result = left_factor.squash_product(right_factor)
        #     self._ensure_valid_rc_graph(result, context="push_factor_through single factor left_squash")
        #     return self._factor_to_normal_tuple(result)

        push_count = 0
        while push_count < len(key) and len(key[push_count]) <= len(left_grass):
            push_count += 1

        # if push_count == 0:
        #     return (left_grass, *key)
        if push_count != 0:

            prefix = key[:push_count]
            suffix = key[push_count:]

            working_rc = self.key_to_rc_graph(prefix)

            if len(working_rc) < len(left_grass):
                working_rc = working_rc.resize(len(left_grass))
                self._ensure_valid_rc_graph(working_rc, context="push_factor_through enlarge prefix rc")

            # Prefix squash, then left-squash by the incoming Grass factor.
            # twisted_rc = working_rc.left_squash(left_grass)
            twisted_rc = left_grass.left_squash(working_rc)
            self._ensure_valid_rc_graph(twisted_rc, context="push_factor_through twisted result")
        else:
            #suffix = key
            #twisted_rc = left_grass
            twisted_rc = left_grass
            for suffix_factor in key:
                if len(suffix_factor) > len(twisted_rc):
                    twisted_rc = twisted_rc.resize(len(suffix_factor))
                    self._ensure_valid_rc_graph(twisted_rc, context="push_factor_through enlarge twisted rc in suffix loop")
                twisted_rc = twisted_rc.squash_product(suffix_factor)
            return self._factor_to_normal_tuple(twisted_rc)

        twisted_rc_base, twisted_rc_grass = twisted_rc.squash_decomp()
        if twisted_rc_base.perm.inv == 0 and twisted_rc_grass.perm.inv == 0:
            return suffix
        if twisted_rc_base.perm.inv == 0:
            return self.push_factor_through(twisted_rc_grass, suffix)
        if twisted_rc_grass.perm.inv == 0:
            return (twisted_rc_base.normalize(), *suffix)
        return (twisted_rc_base.normalize(), *self.push_factor_through(twisted_rc_grass, suffix))

    def _is_strictly_increasing_length_tuple(self, key: tuple) -> bool:
        return all(len(key[i]) < len(key[i + 1]) for i in range(len(key) - 1))

    def _identity_free(self, key: tuple) -> tuple:
        return tuple(rc for rc in key if rc.perm.inv != 0)

    def _factor_to_normal_tuple(self, rc: RCGraph) -> tuple:
        """Maximally factor a single RC graph using repeated squash_decomp on the leftmost factor."""
        self._ensure_valid_rc_graph(rc, context="_factor_to_normal_tuple input")
        if rc.perm.inv == 0:
            return ()
        tup = (rc,)
        if rc.is_full_grass:
            return tup
        if max(rc.perm.descents()) < len(rc) - 1:
            tup = (rc.normalize(),)
        old_tup = tup
        while not tup[0].is_full_grass:
            nepple = tup[0].squash_decomp()
            tup = (nepple[0].normalize(), nepple[1].normalize(), *tup[1:])
            if len(tup) == len(old_tup):
                raise ValueError(f"Squash decomposition did not reduce size: {tup}")
            old_tup = tup

        for factor in tup:
            self._ensure_valid_rc_graph(factor, context="_factor_to_normal_tuple output")
        assert self.key_to_rc_graph(tup) == rc.normalize(), f"Factorization tuple does not evaluate back to original RC: {tup} evaluates to {self.key_to_rc_graph(tup)}, expected {rc.normalize()}"
        return self._identity_free(tup)

    def _left_multiply_factor(self, left_elem: RCGraph, key: tuple) -> tuple:
        """Multiply one Grassmannian factor on the left of a normal-form key."""
        return self._insert_factor_smash(left_elem, key)

    def _insert_factor_smash(self, left_elem: RCGraph, key: tuple) -> tuple:
        """Single-step smash-style insertion of one Grass factor into a normal-form key."""
        self._ensure_valid_rc_graph(left_elem, context="_left_multiply_factor input")
        if left_elem.perm.inv == 0:
            return key

        if not _is_full_grassmannian_rc(left_elem):
            raise ValueError(f"Left factor must be full Grassmannian, got: {left_elem}")

        if len(key) == 0:
            resized = left_elem.resize(len(left_elem))
            self._ensure_valid_rc_graph(resized, context="left multiply empty key resize")
            return self._factor_to_normal_tuple(resized)

        # first = key[0]
        # k = _descent_of_grass(left_elem)
        # p = _descent_of_grass(first)

        # if p >= k:
        #     resized = left_elem.resize(len(left_elem))
        #     self._ensure_valid_rc_graph(resized, context="prepend case resize")
        #     return self._factor_to_normal_tuple((resized, *key))

        # if p == k:
        #     merged = self.same_size_product(left_elem, first)
        #     return (merged, *key[1:])

        return self.push_factor_through(left_elem, key)

    def _normalize_key(self, key: tuple) -> tuple:
        """Normalize an RCGraph tuple to strict increasing-length Grass normal form."""
        key = self._identity_free(key)
        result: tuple = ()
        for factor in reversed(key):
            if not _is_full_grassmannian_rc(factor):
                expanded = self._factor_to_normal_tuple(factor)
                for part in reversed(expanded):
                    result = self._left_multiply_factor(part, result)
            else:
                result = self._left_multiply_factor(factor.resize(len(factor)), result)
        return self._identity_free(result)

    def _mul_keys_legacy(self, left_key: tuple, right_key: tuple) -> tuple:
        raise DeprecationWarning("Legacy product is not guaranteed to produce normal form factors; use product_mode='smash' or 'forced' instead.")
        """Legacy factor-by-factor insertion product (kept for backward compatibility)."""
        result = right_key

        for factor in reversed(left_key):
            result = self._left_multiply_factor(factor, result)
        return result

    def _mul_keys_smash(self, left_key: tuple, right_key: tuple) -> tuple:
        """Unique normal form multiplication via full RC reconstruction.

        1. Evaluate right_key to a single RC via squash-product.
        2. Fold each left Grassmannian factor (largest first = reversed left_key):
                     - len(g) >= len(rc): resize rc up if needed, then left_squash.
                     - len(g) <  len(rc): decompose rc to Grass factors and push g through
                         successively (no direct squash by arbitrary rc).
        3. Canonically factorize the final RC.
        """
        # Evaluate right side to a single RC.
        if right_key:
            right_rc = self.key_to_rc_graph(right_key)
            self._ensure_valid_rc_graph(right_rc, context="_mul_keys_smash right_key evaluation")
        else:
            right_rc = RCGraph()

        # Fold each left Grass factor through, largest first.
        for left_factor in reversed(left_key):
            if left_factor.perm.inv == 0:
                continue
            L = len(left_factor)
            R = len(right_rc)
            if L >= R:
                # Resize right_rc to left_factor's size, then left_squash.
                if R < L:
                    right_rc = right_rc.resize(L)
                    self._ensure_valid_rc_graph(right_rc, context="_mul_keys_smash resize right_rc")
                    right_rc = right_rc.left_squash(left_factor)
                else:
                    right_rc = left_factor.squash_product(right_rc)
                # right_rc = left_factor.left_squash(right_rc)
            else:
                # left_factor is smaller: reconstruct right_rc into Grass factors
                # and push through successively rather than squashing by an arbitrary rc.
                right_key_nf = self._factor_to_normal_tuple(right_rc)
                right_key_nf = self.push_factor_through(left_factor, right_key_nf)
                right_rc = self.key_to_rc_graph(right_key_nf)
            self._ensure_valid_rc_graph(right_rc, context="_mul_keys_smash fold step")

        return self._factor_to_normal_tuple(right_rc)

    def _mul_keys_forced(self, left_key: tuple, right_key: tuple) -> tuple:
        """Forced-normalization product.

        This preserves all factors, then applies one canonical normalization pass.
        Useful for experimentation when local rewrite order is not trusted.
        """
        return self._normalize_key((*left_key, *right_key))

    def _mul_keys(self, left_key: tuple, right_key: tuple) -> tuple:
        if self.product_mode == "smash":
            smash = self._mul_keys_smash(left_key, right_key)
            if self.compare_modes:
                forced = self._mul_keys_forced(left_key, right_key)
                if smash != forced:
                    raise ValueError(f"Smash/forced product mismatch for {left_key} * {right_key}: {smash} != {forced}")
            return smash

        if self.product_mode == "forced":
            forced = self._mul_keys_forced(left_key, right_key)
            if self.compare_modes:
                legacy = self._mul_keys_legacy(left_key, right_key)
                if legacy != forced:
                    raise ValueError(f"Legacy/forced product mismatch for {left_key} * {right_key}: {legacy} != {forced}")
            return forced

        legacy = self._mul_keys_legacy(left_key, right_key)
        if self.compare_modes:
            forced = self._mul_keys_forced(left_key, right_key)
            if legacy != forced:
                raise ValueError(f"Legacy/forced product mismatch for {left_key} * {right_key}: {legacy} != {forced}")
        return legacy

    def _coerce_key(self, key) -> tuple:
        if isinstance(key, RCGraph):
            key = (key,)
        if not isinstance(key, tuple):
            raise TypeError(f"Grass tensor key must be a tuple of RCGraph, got {type(key)}")
        for rc in key:
            if not isinstance(rc, RCGraph):
                raise TypeError(f"Grass tensor factors must be RCGraph, got {type(rc)}")
            self._ensure_valid_rc_graph(rc, context="coerce key")
            if not _is_full_grassmannian_rc(rc):
                raise ValueError(f"Non-full-Grassmannian factor in key: {rc}")
        normalized = self._normalize_key(key)
        if not self._is_strictly_increasing_length_tuple(normalized):
            raise ValueError(f"Key is not in strict increasing-length normal form: {normalized}")
        # if not self._is_fully_factorized_key(normalized):
        #     raise ValueError(f"Key is not fully factorized (length slots/permutation bound): {normalized}")
        return normalized

    def from_dict(self, dct):
        accum = {}
        for key, coeff in dct.items():
            if coeff == S.Zero:
                continue
            nkey = self._coerce_key(key)
            accum[nkey] = accum.get(nkey, S.Zero) + coeff
        return self.dtype({k: v for k, v in accum.items() if v != S.Zero})

    def __call__(self, key):
        return self.from_dict({self._coerce_key(key): S.One})

    def from_rc_graph(self, rc: RCGraph):
        """Create a single-key element from an arbitrary RC graph via normal-form factorization."""
        if not isinstance(rc, RCGraph):
            raise TypeError(f"from_rc_graph expects RCGraph, got {type(rc)}")
        key = self._factor_to_normal_tuple(rc)
        assert rc == self.key_to_rc_graph(key), f"Factorization mismatch: {rc} != {self.key_to_rc_graph(key)}"
        return self.from_dict({key: S.One})

    def mul(self, a, b):
        if not isinstance(b, GrassTensorAlgebraElement):
            return super().mul(a, b)
        result = {}
        for left_key, left_coeff in a.items():
            for right_key, right_coeff in b.items():
                key = self._mul_keys(left_key, right_key)
                result[key] = result.get(key, S.Zero) + left_coeff * right_coeff
        return self.from_dict(result)

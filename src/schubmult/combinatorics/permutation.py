"""The core indexing object used throughout ``schubmult``: finite permutations.

A `Permutation` wraps a 1-indexed array (``self[i]`` is the value at 0-indexed
position ``i``, extended by the identity ``self[i] = i + 1`` past its window)
and is immutable/hashable/cached (equal permutations of different window
lengths, e.g. ``[2, 1]`` and ``[2, 1, 3]``, compare and hash equal).
Multiplication ``p * q`` is composition of the underlying functions,
``(p * q)[i] = p[q[i] - 1]``; ``~p`` is the inverse. Most of the combinatorics
library is built from a `Permutation`'s Lehmer code (``code``/``trimcode``),
inversions, and Bruhat/weak order.
"""

import math
from functools import cache, cached_property

import schubmult.utils.logging as lg
import schubmult.utils.perm_utils as sl
from schubmult.utils._printable import LazyPrintable

logger = lg.get_logger(__name__)

zero = 0
n = 100


class Permutation(LazyPrintable):
    """A finite permutation, stored as a 1-indexed array and extended by the identity.

    Construct from array form, e.g. ``Permutation([2, 1, 3])``, or via
    ``uncode(lehmer_code)``. Supports composition (``*``), inversion (``~``),
    Bruhat order (``<=``/``bruhat_leq``), and iteration over the window
    (``list(perm)``, ``perm[i]`` 0-indexed).
    """

    def apply(self, arr):
        """Permute the elements of ``arr`` by this permutation: return ``[arr[self[i] - 1] for i in range(len(arr))]``."""
        return tuple([arr[self[i] - 1] for i in range(len(arr))])

    @property
    def is_reducible(self):
        """Whether ``self`` splits as a direct sum of two smaller permutations across some fixed point."""
        if self.inv == 0:
            return False
        for i in range(1, len(self)):
            if self[i - 1] == i:
                parabolic_list = list(range(1, i)) + list(range(i + 1, len(self) + 1))
                if self.min_coset_rep(*parabolic_list).inv == 0:
                    return True
        return False

    def reduce(self, start_spot=1, strict=False):
        """Split ``self`` at a fixed point into ``(perm1, perm2)`` on disjoint blocks of values, or ``None`` if not reducible."""
        if self.inv == 0:
            return None
        for i in range(start_spot - 1, len(self)):
            if self[i - 1] == i:
                parabolic_list = list(range(1, i)) + list(range(i + 1, len(self) + 1))
                coset_rep, _residue = self.coset_decomp(*parabolic_list)
                if coset_rep.inv == 0:
                    perm1 = Permutation(self[: i - 1])
                    perm2 = Permutation([a - i for a in self[i:]])
                    return perm1, perm2
            if strict:
                return None
        return None

    def act_root(self, a, b):
        """Image of the root/pair ``(a, b)`` (1-indexed positions) under ``self``: ``(self[a-1], self[b-1])``."""
        return self[a - 1], self[b - 1]

    def one_dominates(self, other):
        """Whether ``self`` one-step dominates ``other`` (one level of the recursive ``dominates`` test)."""
        return _one_dominates(self, other)

    def dominates(self, other):
        """Whether ``self`` dominates ``other`` in the sense used by the dual Pieri / positivity rules."""
        return _dominates(self, other)

    @property
    def antiperm(self):
        """Conjugate of ``self`` by the longest element ``w0``: ``w0 * self * w0``."""
        w0 = Permutation.w0(len(self))
        return w0 * (self) * w0

    def weight_coset_decomp(self, dominant_weight):
        """Coset decomposition of ``self`` with respect to the stabilizer of ``dominant_weight``; see ``coset_decomp``."""
        fixers = Permutation.fixers(dominant_weight)
        return self.coset_decomp(*list(fixers))

    def min_of_weight_coset(self, dominant_weight):
        """Minimal-length coset representative of ``self`` in the stabilizer of ``dominant_weight``."""
        return self.weight_coset_decomp(dominant_weight)[0]

    def max_of_weight_coset(self, dominant_weight):
        """Maximal-length coset representative of ``self`` in the stabilizer of ``dominant_weight``."""
        return self.weight_coset_decomp(dominant_weight)[0] * Permutation.longest_element(*Permutation.fixers(dominant_weight))

    # low_perm is sorting_perm of lowest weight, reverse=True
    @staticmethod
    def does_demazure_crystal_tensor_decompose(dominant_weight1, low_perm1, dominant_weight2, low_perm2):
        """Whether the tensor product of the two Demazure crystals (given by their dominant weights and
        lowest-weight sorting permutations) decomposes as their `does_demazure_crystal_tensor_decompose`
        criterion predicts: the minimal coset representative of ``low_perm1`` in weight ``dominant_weight1``
        lies below the longest element of the left descents of the maximal coset representative of
        ``low_perm2`` in weight ``dominant_weight2``, in Bruhat order.
        """
        # min dom weight 1 low perm 1 is in coset of max dom weight 2 low perm 2
        min_v = low_perm1.min_of_weight_coset(dominant_weight1)
        max_w = low_perm2.max_of_weight_coset(dominant_weight2)

        left_descents = (~max_w).descents(zero_indexed=False)
        return min_v.bruhat_leq(Permutation.longest_element(*left_descents))

    def coset_decomp(self, *descs):
        """Decompose ``self = reduced_perm * w_J`` where ``w_J`` lies in the parabolic subgroup generated
        by the simple reflections at 1-indexed positions ``descs`` and ``reduced_perm`` is the minimal-length
        coset representative (no descents in ``descs``).

        Returns:
            tuple: ``(reduced_perm, w_J)``.
        """
        descs = set(descs)
        reduced_perm = self
        w_J = Permutation([])
        found = True
        while found:
            found = False
            for d in reduced_perm.descents():
                if d + 1 in descs:
                    w_J = ~((~w_J).swap(d, d + 1))
                    reduced_perm = reduced_perm.swap(d, d + 1)
                    found = True
                    break
        return reduced_perm, w_J

    def min_coset_rep(self, *descs):
        """Minimal-length coset representative of ``self`` for the parabolic subgroup generated at ``descs``."""
        return self.coset_decomp(*descs)[0]

    def max_coset_rep(self, *descs):
        """Maximal-length coset representative of ``self`` for the parabolic subgroup generated at ``descs``."""
        red, _w_J = self.coset_decomp(*descs)
        return red * Permutation.longest_element(*descs)

    @classmethod
    def longest_element(cls, *descs):
        """Longest element of the parabolic subgroup generated by the simple reflections at 1-indexed positions ``descs``."""
        perm = Permutation([])
        did_one = True
        while did_one:
            did_one = False
            for i in range(len(descs)):
                j = descs[i] - 1
                if perm[j] < perm[j + 1]:
                    perm = perm.swap(j, j + 1)
                    did_one = True
        return perm

    @classmethod
    def w0(cls, n):
        """The longest element of the symmetric group on ``n`` letters: ``[n, n-1, ..., 1]``."""
        return cls.from_code([n - 1 - i for i in range(n - 1)])

    @classmethod
    @cache
    def all_permutations(cls, n):
        """All permutations of ``{1, ..., n}``, as a list of `Permutation`."""
        from itertools import permutations

        return [cls(perm) for perm in list(permutations(list(range(1, n + 1))))]

    def parabolic_reduce(self, *descs):
        """Decompose ``self = reduced_perm * w_J`` where ``w_J`` is generated by the simple reflections
        *not* in ``descs`` and ``reduced_perm`` has no descents outside ``descs``.

        Returns:
            tuple: ``(reduced_perm, w_J)``.
        """
        descs = set(descs)
        reduced_perm = self
        w_J = Permutation([])
        found = True
        while found:
            found = False
            for d in reduced_perm.descents():
                if d not in descs:
                    w_J = ~((~w_J).swap(d, d + 1))
                    reduced_perm = reduced_perm.swap(d, d + 1)
                    found = True
                    break
        return reduced_perm, w_J

    def __truediv__(self, other):
        """Returns a tuple of (self, other) if other is a permutation. Intended for skew elements."""
        if isinstance(other, Permutation):
            return (self, other)
        raise NotImplementedError("Division by non-permutation is not implemented.")

    def __new__(cls, perm):
        return Permutation.__xnew_cached__(cls, tuple(perm))

    print_as_code = False

    @staticmethod
    def fixers(dominant_weight):
        """1-indexed positions ``i`` where ``dominant_weight[i-1] == dominant_weight[i]``
        (the simple reflections fixing ``dominant_weight``, i.e. generating its stabilizer).
        """
        fixers = set()
        for i in range(len(dominant_weight) - 1):
            if dominant_weight[i] == dominant_weight[i + 1]:
                fixers.add(i + 1)
        return fixers

    @classmethod
    def ref_product(cls, *args):
        """Product of the simple reflections ``s_a`` for ``a`` in ``args``, applied left to right."""
        p = cls([])
        for a in args:
            p = p.swap(a - 1, a)
        return p

    @classmethod
    def hecke_ref_product(cls, *args):
        """Like ``ref_product``, but each ``s_a`` is applied only if it is a Bruhat ascent
        (the 0-Hecke / Demazure product of the simple reflections).
        """
        p = cls([])
        for a in args:
            if p[a - 1] < p[a]:
                p = p.swap(a - 1, a)
        return p

    @property
    def code_word(self):
        """A canonical reduced word for ``self``, read off from ``trimcode``."""
        cd = self.trimcode
        word = []
        for i in range(len(cd)):
            word += list(range(i + cd[i], i, -1))
        return tuple(word)

    @property
    def inverse_code_word(self):
        """A canonical reduced word for ``~self``, read off from ``(~self).trimcode``."""
        cd = (~self).trimcode
        word = []
        for i in range(len(cd)):
            word = list(range(i + 1, i + cd[i] + 1)) + word
        return tuple(word)

    def root_swap(self, root):
        """Multiply ``self`` by the reflection swapping the 1-indexed positions ``root = (a, b)``."""
        return self.swap(root[0] - 1, root[1] - 1)

    @classmethod
    def reflection(cls, root):
        """The transposition swapping the 1-indexed positions ``root = (a, b)``."""
        return cls([]).swap(root[0] - 1, root[1] - 1)

    def right_root_at(self, index, word=None):
        """The positive root sent to a negative root by the ``index``-th letter of ``word``
        (default ``self.code_word``), read from the right (post-multiplied by the remaining suffix).
        """
        if word is None:
            word = [*self.code_word]
        return Permutation._right_root_at(index, word)

    def left_root_at(self, index, word=None):
        """Like ``right_root_at``, but the root is transported by the prefix of ``word`` before ``index``."""
        if word is None:
            word = [*self.code_word]
        return Permutation._left_root_at(index, word)

    @staticmethod
    def _right_root_at(index, word):
        word_piece = word[index + 1 :]
        apply = ~Permutation.ref_product(*word_piece)
        root = apply.act_root(word[index], word[index] + 1)
        return root

    # @staticmethod
    # def _hecke_right_root_at(index, word):
    #     word_piece = word[index + 1 :]
    #     apply = ~Permutation.ref_product(*word_piece)
    #     root = apply.act_root(word[index], word[index] + 1)
    #     return root

    @staticmethod
    def _left_root_at(index, word):
        word_piece = word[:index]
        apply = Permutation.ref_product(*word_piece)
        root = apply.act_root(word[index], word[index] + 1)
        return root

    @cache
    def all_reduced_words(self):
        """All reduced words of `self`, by peeling descents."""
        if self.inv == 0:
            return {()}
        out = set()
        for d in self.descents():
            sub = self.swap(d, d + 1)
            for w in sub.all_reduced_words():
                out.add((*w, d + 1))
        return out

    @staticmethod
    def all_reduced_subwords(word):
        """All reduced subwords of `self`, by peeling descents."""
        if len(word) == 0:
            return {()}
        out = set()
        d = word[-1]
        sub = word[:-1]
        for w in Permutation.all_reduced_subwords(sub):
            out.add(w)
            check_perm = Permutation.ref_product(*w)
            if check_perm[d-1] < check_perm[d]:
                out.add((*w, d))
        return out

    @staticmethod
    def all_subwords(word):
        """All subwords of `self`, by peeling descents."""
        if len(word) == 0:
            return {()}
        out = set()
        d = word[-1]
        sub = word[:-1]
        subwords = Permutation.all_subwords(sub)
        out.update(subwords)
        for w in subwords:
            out.add((*w, d))
        return out

    @staticmethod
    def commutation_class_of(word):
        """All words obtainable from ``word`` by commuting adjacent far-apart letters (``|a - b| >= 2``)."""
        stack = [tuple(word)]
        ret = set()
        while len(stack) > 0:
            u = stack.pop()
            ret.add(u)
            for i in range(len(u) - 1):
                if abs(u[i] - u[i + 1]) >= 2:
                    v = (*u[:i], u[i + 1], u[i], *u[i + 2:])
                    if v not in ret:
                        stack.append(v)
        return ret

    @staticmethod
    def forest_class_of(word):
        """Words reachable from ``word`` by commutation moves that also preserve the indexed-forest
        insertion structure (``omega_insertion``); a refinement of ``commutation_class_of``.
        """
        from .indexed_forests import omega_insertion, word_to_pairinj_labeled

        stack = [tuple(word)]
        ret = set()
        while len(stack) > 0:
            u = stack.pop()
            ret.add(u)
            for i in range(len(u) - 1):
                if abs(u[i] - u[i + 1]) >= 2:
                    da_word = word_to_pairinj_labeled(tuple(reversed(u[i:])))
                    P_right = omega_insertion(da_word[:-2])[0]
                    if len(P_right.separators(da_word[-2], da_word[-1])) > 1:
                        v = (*u[:i], u[i + 1], u[i], *u[i + 2:])
                        if v not in ret:
                            stack.append(v)
        return ret

    def code_index_of_index(self, index):
        """The position in ``trimcode`` whose block of ``code_word`` letters contains position ``index``."""
        running_sum = 0
        running_code_index = 0
        for code_index, code_elem in enumerate(self.trimcode):
            if code_elem == 0:
                continue
            running_sum += code_elem
            if running_sum > index:
                return running_code_index
            running_code_index += 1
        return len(self.trimcode)

    @staticmethod
    def cycle(p, q):
        """
        Construct the cycle permutation used elsewhere in the code.
        Kept as a staticmethod on Permutation for call sites like Permutation.cycle(p,q).
        """
        return Permutation(list(range(1, p)) + [i + 1 for i in range(p, p + q)] + [p])

    @staticmethod
    @cache
    def __xnew_cached__(_class, perm):
        return Permutation.__xnew__(_class, perm)

    @staticmethod
    def __xnew__(_class, perm):
        p = tuple([int(x) for x in sl.permtrim_list([*perm])])
        # s_perm = spp.Permutation([i - 1 for i in p])
        obj = object.__new__(_class)
        obj._args = (p,)
        # obj._s_perm = tuple([i - 1 for i in p])
        obj._perm = p
        obj._hash_code = hash(p)
        cd = sl.old_code(p)
        obj._unique_key = (len(p), sum([cd[i] * math.factorial(len(p) - 1 - i) for i in range(len(cd))]))
        return obj

    @cached_property
    def _arr(self):
        """Cached numpy array representation (1-indexed values)."""
        import numpy as np

        return np.array(self._perm, dtype=int)

    @property
    def args(self):
        return self._args

    @classmethod
    def sorting_perm(cls, itera, reverse=False):
        """The permutation that sorts ``itera`` into (by default) increasing order."""
        L = [i + 1 for i in range(len(itera))]
        L.sort(key=lambda i: itera[i - 1], reverse=reverse)
        return Permutation(L)

    def right_act(self, lst):
        """Permute the entries of ``lst`` by this permutation, preserving ``lst``'s type (list stays a list)."""
        if isinstance(lst, list):
            return [lst[self[i] - 1] for i in range(len(lst))]
        return tuple([lst[self[i] - 1] for i in range(len(lst))])

    def bruhat_leq(perm, perm2):
        """Whether ``perm <= perm2`` in Bruhat order (equivalently, ``perm``'s tableau criterion
        against ``perm2`` on every prefix of their windows).
        """
        if perm.inv == perm2.inv:
            return perm == perm2
        if perm.inv > perm2.inv:
            return False
        ml = max(len(perm), len(perm2))
        full_perm = [perm[i] for i in range(ml)]
        full_perm2 = [perm2[i] for i in range(ml)]
        for i in range(1, ml):
            arr1 = list(full_perm[:i])
            arr2 = list(full_perm2[:i])
            arr1.sort()
            arr2.sort()
            if any(a1 > a2 for a1, a2 in zip(arr1, arr2)):
                return False
        return True

    @classmethod
    def from_code(cls, cd):
        """Alias for ``uncode(cd)``: the permutation with Lehmer code ``cd``."""
        return uncode(cd)

    # def _latex(self, printer):
    #     if Permutation.print_as_code:
    #         return printer._print(self.trimcode)
    #     return printer._print(list(self._perm))

    # pattern is a list, not a permutation
    def has_pattern(self, pattern):
        """Whether ``self`` contains ``pattern`` (a plain list, not necessarily reduced) as a pattern:
        some subsequence of ``self``'s window order-isomorphic to ``pattern``.
        """
        if self == Permutation(pattern):
            return True
        if len(self._perm) <= len(Permutation(pattern)):
            return False
        expanded = list(self) + list(range(len(self) + 1, len(pattern) + 1))
        for i in range(len(expanded)):
            rmval = expanded[i]
            perm2 = [*expanded[:i], *expanded[i + 1 :]]
            perm2 = tuple([val - 1 if val > rmval else val for val in perm2])
            if Permutation(perm2).has_pattern(pattern):
                return True
        return False

    def _pretty(self, printer=None):
        return printer._print_Tuple(tuple(self))

    def _sympystr(self, printer=None):
        from sympy.printing.str import StrPrinter

        if printer is None:
            printer = StrPrinter()
        if Permutation.print_as_code:
            return printer.doprint(self.trimcode)
        return printer.doprint(tuple(self._perm))

    def _latex(self, printer):
        if Permutation.print_as_code:
            return printer.doprint(self.trimcode)
        return printer.doprint(tuple(self._perm))

    def __call__(self, *tup):
        if len(tup) == 1:
            if isinstance(tup[0], list | tuple):
                tup = tup[0]
            else:
                return self._perm[tup[0] - 1]
        return tuple(self[i - 1] for i in tup)

    def zero_indexed_descents(self):
        """0-indexed descent positions: ``i`` such that ``self[i] > self[i+1]``."""
        desc = set()
        for i in range(len(self._perm) - 1):
            if self[i] > self[i + 1]:
                desc.add(i)
        return desc

    def descents(self, zero_indexed=True):
        """Descent positions of ``self``, 0-indexed by default or 1-indexed if ``zero_indexed=False``."""
        if zero_indexed:
            return self.zero_indexed_descents()
        return {i + 1 for i in self.zero_indexed_descents()}

    def get_cycles(self, sort_min=False):
        """Cycle decomposition of ``self`` as a list of tuples; ``sort_min`` rotates each cycle to start
        at its minimum element instead of the sympy convention.
        """
        return self.get_cycles_cached(sort_min)

    @cache
    def get_cycles_cached(self, sort_min):
        import sympy.combinatorics.permutations as spp

        if not sort_min:
            return [tuple(sl.cyclic_sort([i + 1 for i in c])) for c in spp.Permutation([k - 1 for k in self._perm]).cyclic_form]
        return [tuple(sl.cyclic_sort_min([i + 1 for i in c])) for c in spp.Permutation([k - 1 for k in self._perm]).cyclic_form]

    @classmethod
    def from_cycles(cls, cycle_iter):
        """Build a permutation from a cycle decomposition (sequence of cycles, each a sequence of 1-indexed values)."""
        import sympy.combinatorics.permutations as spp

        spoing = spp.Permutation(*cycle_iter)
        return cls([a + 1 for a in spoing.array_form])

    @property
    def code(self):
        """Lehmer code of ``self``: ``code[i]`` counts ``j > i`` with ``self[i] > self[j]``."""
        return [*self._cached_code()]

    @cache
    def _cached_code(self):
        return sl.old_code(self._perm)

    @property
    def graph(self):
        """The permutation matrix support as a set of 1-indexed pairs ``{(i, self[i])}``."""
        return {(i + 1, self[i]) for i in range(len(self._perm))}

    def reduced_with(self, other):
        """Whether ``self * other`` is length-additive: ``inv(self) + inv(other) == inv(self * other)``."""
        return (self.inv + other.inv) == (self * other).inv

    @property
    def shape(self):
        """The Lehmer code sorted into weakly decreasing order (a partition)."""
        return tuple(sorted(self.code, reverse=True))

    @cached_property
    def inversion_set(self):
        """Set of inversions ``(i, j)`` (1-indexed, ``i < j``) with ``self[i-1] > self[j-1]``."""
        inv_set = set()
        for i in range(len(self._perm)):
            for j in range(i + 1, len(self._perm)):
                if self[i] > self[j]:
                    inv_set.add((i + 1, j + 1))
        return inv_set

    # left weak order
    def weak_order_leq(self, other):
        """Whether ``self <= other`` in left weak order (``self``'s inversion set is a subset of ``other``'s)."""
        return self.inversion_set.issubset(other.inversion_set)

    def weak_order_meet(self, other):
        """Meet (greatest lower bound) of ``self`` and ``other`` in left weak order."""
        if self.weak_order_leq(other):
            return self
        if other.weak_order_leq(self):
            return other
        invset = self.inversion_set.intersection(other.inversion_set)
        descents = {a for (a, b) in invset if a + 1 == b}
        if len(descents) == 0:
            return Permutation([])
        a = max(descents)
        downself = self.swap(a - 1, a)
        downother = other.swap(a - 1, a)
        return downself.weak_order_meet(downother).swap(a - 1, a)

    def weak_order_join(self, other):
        """Join (least upper bound) of ``self`` and ``other`` in left weak order, via ``w0``-duality with the meet."""
        max_len = max(len(self), len(other))
        w0 = Permutation.w0(max_len)
        return (self * w0).weak_order_meet(other * w0) * w0

    @property
    def diagram(self):
        """The Rothe diagram of ``self``: cells ``(i, j)`` with ``self[i-1] > j`` and ``(~self)[j-1] > i``."""
        diag = set()
        for i in range(len(self._perm)):
            for j in range(len(self._perm)):
                if self[i] > j + 1 and (~self)[j] > i + 1:
                    diag.add((i + 1, j + 1))
        return diag

    @property
    def rothe_diagram(self):
        """The graph of ``self`` as a set of 1-indexed pairs (see ``diagram`` for the Rothe diagram cells)."""
        return {(i + 1, self[i]) for i in range(len(self))}

    @cached_property
    def max_descent(self):
        """Number of entries in ``trimcode`` (one past the last nonzero code entry)."""
        return len(self.trimcode)

    @property
    def maximal_corner(self):
        """The maximal corner ``(maxd, end_spot)`` of ``self``'s diagram, used by ``pivots``/``pivot_transition``."""
        maxd = len(self.trimcode)
        end_spot = max(self[i] for i in range(maxd, len(self)) if self[i] < self[maxd - 1])
        return (maxd, end_spot)

    @classmethod
    def from_partial(cls, partial_perm):
        """Complete a partial assignment (a list with some ``None`` entries) to a full permutation,
        filling the gaps with the missing values in increasing order.
        """
        max_required = max([a for a in partial_perm if a is not None], default=len(partial_perm))
        partial_perm = list(partial_perm)
        if len(partial_perm) < max_required:
            partial_perm += [None] * (max_required - len(partial_perm))
        search_space = {i for i in partial_perm if i is not None}

        # Need enough values to fill all None positions
        full_perm = [i + 1 for i in range(max(max_required, len(partial_perm))) if i + 1 not in search_space]
        perm = [*partial_perm]
        j = 0
        for i in range(len(perm)):
            if perm[i] is None:
                perm[i] = full_perm[j]
                j += 1
        return cls(perm)

    @cache
    def pivots(self, a=None, b=None):
        """Return the set of pivot positions for a maximal corner (a,b)."""
        if a is None or b is None:
            a, b = self.maximal_corner
        piv = set()
        for i in range(1, len(self) + 1):
            j = self[i - 1]
            if i >= a or j >= b:
                continue
            good = True
            for i_prime in range(i, a + 1):
                if not good:
                    break
                for j_prime in range(j, b + 1):
                    if (i, j) == (i_prime, j_prime) or (i_prime, j_prime) == (a, b):
                        continue
                    if self[i_prime - 1] == j_prime:
                        good = False
                        break
            if good:
                piv.add(i)
        return piv

    def pivot_transition(self, pivot_set):
        """Grothendieck transition for a given pivot set at the maximal corner. Returns the resulting permutation."""
        if self.inv == 0:
            raise ValueError("Cannot perform pivot transition on the identity permutation.")
        if not pivot_set.issubset(self.pivots()):
            raise ValueError(f"Invalid pivot set {pivot_set} for permutation {self} with pivots {self.pivots()}")
        pivot_list = sorted(pivot_set, reverse=True)
        maxd, b = self.maximal_corner
        cycle = [maxd, *pivot_list]
        if len(pivot_list) == 0:
            cycle_perm = Permutation([])
        else:
            cycle_arr = list(range(1, len(self) + 1))
            for i in range(len(cycle) - 1):
                cycle_arr[cycle[i] - 1] = cycle[i + 1]
            cycle_arr[cycle[-1] - 1] = cycle[0]
            cycle_perm = Permutation(cycle_arr)
        b_prime = (~self)[b - 1]
        return self.swap(maxd - 1, b_prime - 1) * cycle_perm

    def pad_code(self, length):
        """``trimcode`` padded with trailing zeros to ``length``."""
        if length < len(self.trimcode):
            raise ValueError("Cannot pad to a length shorter than the trimcode")
        return tuple(list(self.trimcode) + [0 for i in range(length - len(self.trimcode))])

    @cached_property
    def trimcode(self):
        """Lehmer code truncated to drop trailing zeros (length equals the last descent position)."""
        if self._perm == ():
            return []
        return self.code[: max(self.descents(False), default=0)]

    def mul_dominant(self):
        """Left factor of ``self`` in its decomposition against the minimal dominant permutation above it."""
        return ~((~self).minimal_dominant_above())

    def strict_mul_dominant(self, size=None):
        """Variant of ``mul_dominant`` built from the strict theta (strictly decreasing dominant code)."""
        if size is None:
            return uncode((~(uncode((~self).theta()))).strict_theta())
        the_perm = uncode([self.trimcode[a] + 1 if a < len(self.trimcode) else 1 for a in range(size)])
        return uncode((~(uncode((~the_perm).theta()))).strict_theta())

    def shiftup(self, k):
        """``self`` with its Lehmer code shifted right by ``k`` (``k`` leading zero code entries prepended)."""
        return Permutation.from_code(k * [0] + self.code)

    @cached_property
    def inv(self):
        """Length of ``self``: the number of inversions, i.e. ``sum(self.code)``."""
        return sum(self.code)

    @property
    def is_dominant(self):
        """Whether ``self`` equals the minimal dominant permutation above it (its code is already weakly decreasing)."""
        return self.minimal_dominant_above() == self

    @property
    def is_strict_dominant(self):
        """Whether ``self`` is dominant with a strictly decreasing ``trimcode``."""
        return self.is_dominant and all(self.trimcode[i] > self.trimcode[i + 1] for i in range(self.max_descent - 1))

    @property
    def is_vexillary(self):
        """Whether ``self`` avoids the pattern ``2143`` (equivalently, its Schubert polynomial is a
        single Schur polynomial in the ``trimcode``-shape).
        """
        return not self.has_pattern([2, 1, 4, 3])

    def __reduce__(self):
        return (self.__class__, (self._perm,))

    def swap(self, i, j):
        """Multiply ``self`` by the transposition of 0-indexed positions ``i`` and ``j`` (window extended as needed)."""
        if i > j:
            i, j = j, i
        if j >= len(self._perm):
            # Need to extend - fall back to list operations
            new_perm = [*self._perm]
            new_perm.extend(range(len(new_perm) + 1, j + 2))
            new_perm[i], new_perm[j] = new_perm[j], new_perm[i]
            return Permutation(new_perm)
        # Fast path using numpy
        new_arr = self._arr.copy()
        new_arr[i], new_arr[j] = new_arr[j], new_arr[i]
        return Permutation(new_arr)

    def rslice(self, start, stop):
        """Window values at 0-indexed positions ``start`` (inclusive) to ``stop`` (exclusive), extended by fixed points."""
        ttup = [*self._perm, *list(range(len(self._perm) + 1, stop + 2))]
        return ttup[start:stop]

    def __getitem__(self, i):
        try:
            return self._perm[i]
        except Exception:
            if isinstance(i, slice):
                return [self[ii] for ii in range(i.start if i.start is not None else 0, i.stop if i.stop is not None else len(self))]
            if i >= len(self._perm):
                return i + 1

    def __setitem__(self, i, v):
        raise NotImplementedError

    def __hash__(self):
        return self._hash_code

    def __matmul__(self, other):
        """Demazure product"""
        word = other.code_word
        ret = self
        if len(word) == 0:
            return ret
        for letter in word:
            if ret[letter - 1] < ret[letter]:
                ret = ret.swap(letter - 1, letter)
        return ret

    def __mul__(self, other):
        a, b = self._perm, other._perm
        la = len(a)
        if len(b) < la:
            b = b + tuple(range(len(b) + 1, la + 1))
        return Permutation([a[j - 1] if j <= la else j for j in b])

    def __iter__(self):
        yield from self._perm.__iter__()

    def __getslice__(self, i, j):
        return self._perm[i:j]

    # def __str__(self):
    #     return str(self._perm)

    def __add__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*self._perm, *other]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    # def _sympyrepr(self, printer):
    #     return f"Permutation({list(self._perm)})"

    def __radd__(self, other):
        if not isinstance(other, list):
            raise NotImplementedError
        permlist = [*other, *self._perm]
        try:
            return Permutation(permlist)
        except Exception:
            return permlist

    def __eq__(self, other):
        if isinstance(other, Permutation):
            # print(f"{other._perm= } {self._perm=} {type(self._perm)=}")
            # return other._perm == self._perm
            return other._unique_key == self._unique_key
        if isinstance(other, list):
            # print(f"{[*self._perm]= } {other=}")
            return [*self._perm] == other
        if isinstance(other, tuple):
            # print(f"{self._perm=} {other=}")
            return self._perm == other
        return False

    def __len__(self):
        # print("REMOVE THIS")
        return max(len(self._perm), 2)

    def __le__(self, other):
        return self.bruhat_leq(other)

    # def __lt__(self, other):
    #     return self != other and self.bruhat_leq(other)

    def __invert__(self):
        new_arr = [0] * len(self._perm)
        for i, v in enumerate(self._perm, 1):
            new_arr[v - 1] = i
        return Permutation(new_arr)

    def __str__(self):
        # Same output as sstr (StrPrinter prints int tuples/lists like Python) without importing sympy
        if Permutation.print_as_code:
            return f"[{', '.join(map(str, self.trimcode))}]"
        if len(self._perm) == 1:
            return f"({self._perm[0]},)"
        return f"({', '.join(map(str, self._perm))})"

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return tuple(self) < tuple(other)

    def pattern_at(self, *indices):
        """The (inverse sorting) pattern induced by ``self`` on the given 1-indexed ``indices``."""
        indices = sorted(indices)
        seq = [self[i] for i in indices]
        return ~Permutation.sorting_perm(seq)

    def minimal_dominant_above(self):
        """The minimal dominant permutation ``>= self`` in Bruhat order: ``uncode(self.theta())``."""
        return uncode(self.theta())

    @property
    def foundational_root(self):
        """The pivot pair ``(k, mx)`` marking ``self``'s last descent block and the last position dropping below it."""
        if self.inv == 0:
            return None
        mx = -1
        k = max(self.descents()) + 1
        for i in range(k + 1, len(self) + 1):
            if self[k - 1] > self[i - 1]:
                mx = i
        return (k, mx)

    @cache
    def _cached_strict_theta(self):
        ret = [*self.trimcode]
        did_one = True
        while did_one:
            did_one = False
            for i in range(len(ret) - 2, -1, -1):
                if ret[i + 1] != 0 and ret[i] <= ret[i + 1]:
                    ret[i], ret[i + 1] = ret[i + 1] + 1, ret[i]
                    did_one = True
                    break
        while len(ret) > 0 and ret[-1] == 0:
            ret.pop()
        return tuple(ret)

    def theta(self):
        """Dominant (weakly decreasing) sequence bounding ``self``'s code, used throughout the v-path
        multiplication algorithms (see ``schubmult.mult``).
        """
        return [*self._cached_theta()]

    def medium_theta(self):
        """Variant of ``theta`` used by the \"fast\"/merged-layer multiplication kernels."""
        return [*self._cached_medium_theta()]

    def strict_theta(self):
        """Variant of ``theta`` with strictly decreasing entries (no repeated nonzero layers)."""
        return [*self._cached_strict_theta()]

    def maximal_sortable_below(self):
        """The maximal \"sortable\" permutation ``<= self`` (code has no gap of more than 1 between
        consecutive entries), used by ``mul_sortable``.
        """
        return self._cached_maximal_sortable_below()

    def mul_sortable(self):
        """Left factor of ``self`` against its ``maximal_sortable_below`` decomposition (dual to ``mul_dominant``)."""
        return ~((~self).maximal_sortable_below())

    @cache
    def _cached_maximal_sortable_below(self):
        working_perm = self
        while True:
            L = len(working_perm.trimcode)
            cd = [*working_perm.trimcode, 0]
            loc = max([i for i in range(L) if cd[i] - cd[i + 1] > 1], default=-1)
            if loc == -1:
                return working_perm
            working_perm = working_perm.swap(loc, loc + 1)
        raise ValueError("Should not reach here")

    @cache
    def _cached_theta(self):
        cd = list(self.code)
        for i in range(len(cd) - 1, 0, -1):
            for j in range(i - 1, -1, -1):
                if cd[j] < cd[i]:
                    cd[i] += 1
        cd.sort(reverse=True)
        return tuple(cd)

    @cache
    def _cached_medium_theta(self):
        cd = list(self.code)
        found_one = True
        while found_one:
            found_one = False
            for i in range(len(cd) - 1):
                if cd[i] < cd[i + 1]:
                    found_one = True
                    cd[i], cd[i + 1] = cd[i + 1] + 1, cd[i]
                    break
                if cd[i] == cd[i + 1] and cd[i] != 0 and i > 0 and cd[i - 1] <= cd[i] + 1:
                    cd[i] += 1
                    found_one = True
                    break
        return tuple(cd)


def uncode(cd):
    """The permutation whose Lehmer code is ``cd`` (a list of nonnegative integers)."""
    cd2 = [*cd]
    if cd2 == []:
        return Permutation([])
    max_required = max([cd2[i] + i for i in range(len(cd2))])
    cd2 += [0 for i in range(len(cd2), max_required)]
    fullperm = [i + 1 for i in range(len(cd2) + 1)]
    perm = []
    for i in range(len(cd2)):
        perm += [fullperm.pop(cd2[i])]
    perm += [fullperm[0]]
    return Permutation(perm)


def permtrim(perm):
    """Normalize ``perm`` (a plain array) into a `Permutation` (trims trailing fixed points)."""
    return Permutation(perm)


def cycle(p, q):
    # keep a thin module-level wrapper for backwards compatibility
    return Permutation.cycle(p, q)


def phi1(u):
    """Drop the first entry of ``(~u).code`` and re-invert: one step of the ``dominates`` recursion."""
    c_star = (~u).code
    c_star.pop(0)
    # print(f"{uncode(c_star)=}")
    return ~(uncode(c_star))


def split_perms(perms):
    """Split each permutation in ``perms`` (after the first) into two smaller ones across a
    reducible point, whenever a valid split point exists; used to normalize a chain of
    dominant permutations into minimal reducible pieces.
    """
    perms2 = [perms[0]]
    for perm in perms[1:]:
        cd = perm.code
        index = -1
        not_zero = False
        did = False
        for i in range(len(cd)):
            if cd[i] != 0:
                not_zero = True
            elif not_zero and cd[i] == 0:
                not_zero = False
                index = i
                num_zeros_to_miss = 0
                for j in range(index):
                    if cd[j] != 0:
                        num_zeros_to_miss = max(num_zeros_to_miss, cd[j] - (index - 1 - j))
                num_zeros = 0
                for j in range(index, len(cd)):
                    if cd[j] != 0:
                        break
                    num_zeros += 1
                if num_zeros >= num_zeros_to_miss:
                    cd1 = cd[:index]
                    cd2 = [0 for i in range(index)] + cd[index:]
                    perms2 += [
                        uncode(cd1),
                        uncode(cd2),
                    ]
                    did = True
                    break
        if not did:
            perms2 += [perm]
    return perms2


def _one_dominates(u, w):
    c_star_u = (~u).code
    c_star_w = (~w).code

    a = c_star_u[0]
    b = c_star_w[0]

    for i in range(a, b):
        if i >= len(u) - 1:
            return True
        if u[i] > u[i + 1]:
            return False
    return True


def _dominates(u, w):
    u2 = u
    w2 = w
    while u2.inv > 0 and _one_dominates(u2, w2):
        u2 = phi1(u2)
        w2 = phi1(w2)
    if u2.inv == 0:
        return True
    return False


bad_classical_patterns = [Permutation([1, 4, 2, 3]), Permutation([1, 4, 3, 2]), Permutation([4, 1, 3, 2]), Permutation([3, 1, 4, 2])]

ID_PERM = Permutation([])


@cache
def s(i):
    """The simple reflection swapping 1-indexed positions ``i`` and ``i + 1``."""
    return Permutation([*list(range(1, i)), i + 1, i])

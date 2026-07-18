from __future__ import annotations

import itertools
from collections.abc import Sequence
from functools import cache, cached_property
from typing import ClassVar

from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.schubert_monomial_graph import SchubertMonomialGraph
from schubmult.symbolic import Expr, S, prod
from schubmult.utils._grid_print import GridPrint

from .crystal_graph import CrystalGraph

# from schubmult.utils.perm_utils import _is_compatible

@cache
def pivot_transition(perm2, target_d=None):
    build_groth = set()
    if target_d is None:
        d = perm2.max_descent
    else:
        d = target_d
    if perm2.max_descent < d:
        return perm2
    pivots = perm2.pivots()
    for r in range(1, len(pivots) + 1):
        for pivot_set in itertools.combinations(pivots, r):
            pivot_set = set(pivot_set)
            ptrans = perm2.pivot_transition(pivot_set)
            # print("Debug: ", perm2, pivot_set, ptrans)
            build_groth.update(pivot_transition(ptrans, target_d=d))
    return build_groth


def _is_row_root(row: int, root: tuple[int, int]) -> bool:
    return row is None or (root[0] <= row and root[1] > row)


def _is_compatible(compat_seq, word):
    """Check if compat_seq is compatible with word, i.e. for each prefix of word, the corresponding prefix of compat_seq has all distinct entries."""
    if len(compat_seq) != len(word):
        return False
    if len(compat_seq) == 0:
        return True
    if compat_seq[0] > word[0]:
        return False
    for i in range(1, len(word)):
        if compat_seq[i - 1] > compat_seq[i]:
            return False
        if word[i - 1] <= word[i] and compat_seq[i - 1] == compat_seq[i]:
            return False
        if compat_seq[i] > word[i]:
            return False
    return True


class WCGraph(SchubertMonomialGraph, CrystalGraph, GridPrint, tuple):
    """Word-compatible graph.

    Internal representation matches RCGraph: a tuple of rows where row i (0-indexed)
    is a strictly decreasing tuple of positive integers >= i + 1.

    For WCGraph, the associated permutation is the Demazure product of the
    concatenated row word (1-indexed simple reflections).
    """

    _double = False # condition for double raising/lowering operators

    @property
    def args(self) -> tuple:
        return ()

    def _sympyrepr(self, printer=None) -> str:
        rows = list(self)
        if printer is None:
            return f"WCGraph({rows!r})"
        return f"WCGraph({printer._print(rows)})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, WCGraph):
            return NotImplemented
        return tuple(self) == tuple(other)

    def __hash__(self) -> int:
        return hash((tuple(self), "Tweezers"))

    def trans_co_pipe(self):
        # return self._rebuild([tuple(reversed([])) for i, row in enumerate(reversed(self))])
        new_wc = WCGraph([]).resize(2 * len(self))
        for i in range(1, self.rows + 1):
            for j in range(1, self.cols + 1):
                if not self.has_element(i, j):
                    new_wc = new_wc.toggle_ref_at(i + j, j)
        return new_wc

    def __new__(cls, *args):
        if len(args) == 1:
            try:
                if len(args[0]) == 0:
                    args = ()
            except TypeError:
                pass
        new_args = tuple(tuple(arg) for arg in args)
        return WCGraph.__xnew_cached__(cls, *new_args)

    @staticmethod
    @cache
    def __xnew_cached__(_class, *args):
        return WCGraph.__xnew__(_class, *args)

    @staticmethod
    def __xnew__(_class, *args):
        return tuple.__new__(_class, *args)

    def __init__(self, *args: object) -> None:
        pass

    def _rebuild(self, rows=()) -> WCGraph:
        return type(self)(rows)

    @cached_property
    def perm_word(self) -> tuple[int, ...]:
        ret = []
        for row in self:
            ret = [*ret, *row]
        return tuple(ret)

    @cached_property
    def perm(self) -> Permutation:
        perm = Permutation([])
        for row in self:
            for p in row:
                perm = perm @ Permutation.ref_product(p)
        return perm

    @property
    def hecke_perm(self) -> Permutation:
        return self.perm

    @property
    def is_rc(self) -> bool:
        for i, row in enumerate(self):
            for a in row:
                if a < i + 1:
                    return False
        return True

    @property
    def is_reduced(self):
        return len(self.perm_word) == self.perm.inv

    @property
    def is_valid(self) -> bool:
        for i, row in enumerate(self):
            if any(a < i + 1 for a in row):
                return False
            if any(row[j] <= row[j + 1] for j in range(len(row) - 1)):
                return False
        if not _is_compatible(self.compatible_sequence, self.perm_word):
            return False
        return True

    def shiftup(self, shift: int = 1, check_valid=True) -> WCGraph:
        ret = self._rebuild([tuple(a + shift for a in row) for row in self])
        if check_valid:
            assert ret.is_valid
        return ret

    def normalize(self) -> WCGraph:
        return self.resize(self.perm.max_descent)

    def resize(self, new_length: int) -> WCGraph:
        if new_length < len(self):
            return self.rowrange(0, new_length)
        return self.extend(new_length - len(self))

    def rowrange(self, start: int, end: int | None = None) -> WCGraph:
        if end is None:
            end = len(self)
        if start == end:
            return self._rebuild(())
        return self._rebuild([tuple(a - start for a in row) for row in self[start:end]])

    def extend(self, extra_rows: int) -> WCGraph:
        return self._rebuild([*self, *tuple([()] * extra_rows)])

    def toggle_ref_at(self, i: int, j: int) -> WCGraph:
        if i <= 0 or j <= 0:
            raise IndexError()
        new_row = [*self[i - 1]]
        if i + j - 1 in new_row:
            index = new_row.index(i + j - 1)
            new_row = [*new_row[:index], *new_row[index + 1 :]]
        else:
            index = 0
            if len(new_row) > 0:
                while index < len(new_row) and new_row[index] > i + j - 1:
                    index += 1
            new_row.insert(index, i + j - 1)
        return self._rebuild([*self[: i - 1], tuple(new_row), *self[i:]])

    @cache
    def has_element(self, i: int, j: int) -> bool:
        return i <= len(self) and i + j - 1 in self[i - 1]

    @cached_property
    def length_vector(self) -> tuple[int, ...]:
        return tuple(len(row) for row in self)

    @cached_property
    def weight(self) -> tuple[int, ...]:
        wt = []
        for i, row in enumerate(self):
            wt.extend([i + 1] * len(row))
        return tuple(wt)

    @property
    def rows(self) -> int:
        return len(self)

    @property
    def cols(self) -> int:
        return len(self.perm) - 1

    @property
    def width(self) -> int:
        return self.cols

    @property
    def height(self) -> int:
        return self.rows

    @cached_property
    def compatible_sequence(self) -> tuple[int, ...]:
        seq = []
        for i in range(len(self)):
            for _ in range(len(self[i])):
                seq.append(i + 1)
        return tuple(seq)

    @cache
    def to_mbpd(self, n: int | None = None):
        r"""The marked bumpless pipedream ``Psi(RCP(self))`` (paper
        ``writing/mbpd.solve.tex``, Theorem "T: main").

        This composes the trivial ``WCGraph -> RCP`` repackaging with the
        row-unpop bijection ``Psi``.  ``n`` is the ambient grid size (defaults
        to ``len(self.perm) - 1`` padded to fit the graph).

        Cached: WCGraphs are immutable and hashable, so the (self, n) round
        trip is memoized."""
        from schubmult.combinatorics.mbpd import RCP

        if n is None:
            n = max(len(self.perm) - 1, len(self) + 1, max(self.perm_word, default=0) + 1)
        return RCP.from_wcgraph(self, n=n).psi()

    @classmethod
    @cache
    def from_mbpd(cls, mbpd) -> WCGraph:
        r"""Inverse of :meth:`to_mbpd`: the WCGraph ``RCP(Phi(mbpd))`` obtained
        from the row-pop bijection ``Phi`` (paper Theorem "T: main").

        Cached on the (hashable) ``mbpd`` argument."""
        return mbpd.phi().to_wcgraph()

    @property
    def is_principal(self) -> bool:
        from schubmult.combinatorics.permutation import uncode
        return self.perm == uncode(self.length_vector)

    @classmethod
    def grove_wcs(cls, comp, length=None, base_rc=None) -> set[WCGraph]:
        """Grove polynomial of ``comp`` built by inverting the omega insertion.

        We enumerate the compatible set-valued labelings ``kappa`` of the indexed
        forest ``F = weak_composition_to_indfor(comp)`` (the grove definition), and
        realize each labeling as a WCGraph by writing down its (word, compatible
        sequence) pair *explicitly* -- the published inverse of the set-valued omega
        insertion -- and feeding it to ``WCGraph.from_word_compatible``.

        Concretely:

        * The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
          labeling ``Q`` of the principal reduced RC graph form the omega pair for
          ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
          reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
          so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
          word position ``len(W) - Q(v)`` and carries the reduced letter
          ``ell(v) = W[len(W) - Q(v)]``.
        * A labeling ``kappa`` places the letter ``ell(v)`` into every row
          ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
          ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
          a word whose rows are strictly decreasing, i.e. exactly the data
          ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
          values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
          equals ``x^kappa``.

        Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
        the degree ``-1`` homogenizer recording extra labels beyond one per node.
        """
        from itertools import combinations

        from schubmult.combinatorics.indexed_forests import omega_reduced_word_from_labelings, weak_composition_to_indfor
        from schubmult.combinatorics.permutation import uncode
        from schubmult.combinatorics.rc_graph import RCGraph

        if length is None:
            length = len(comp)
        forest = weak_composition_to_indfor(comp)

        def labelings_below(node, min_value):
            results = []
            for subset_size in range(1, node.rho - min_value + 2):
                for subset in combinations(range(min_value, node.rho + 1), subset_size):
                    top = subset[-1]
                    left_opts = labelings_below(node.left, top) if node.left is not None else [{}]
                    right_opts = labelings_below(node.right, top + 1) if node.right is not None else [{}]
                    for left in left_opts:
                        for right in right_opts:
                            labeling = {node.index: tuple(subset)}
                            labeling.update(left)
                            labeling.update(right)
                            results.append(labeling)
            return results

        all_labelings = [{}]
        for root in forest._roots:
            root_labelings = labelings_below(root, 1)
            all_labelings = [{**base, **choice} for base in all_labelings for choice in root_labelings]

        # Omega pair (P, Q) for F, read off the principal reduced RC graph.
        if base_rc is not None:
            principal = base_rc.resize(length)
        else:
            principal = RCGraph.principal_rc(uncode(comp), length)
        p_labeling, q_labeling = principal.omega_invariant

        # Gamma (the published inverse of the insertion) reconstructs the reduced word
        # from (P, Q); reversing matches the forward convention of perm_word.
        word = list(reversed(omega_reduced_word_from_labelings(p_labeling, q_labeling)))

        # The reduced letter carried by node v is W at v's Q-position.
        letter_of = {node.index: word[len(word) - q_labeling(node.index)] for node in p_labeling.forest.inorder_traversal}

        wc_set = set()
        for labeling in all_labelings:
            # Explicit (word, compatible sequence): node v's letter into each row of kappa(v).
            entries = [(row, letter_of[index]) for index, label_set in labeling.items() for row in label_set]
            entries.sort(key=lambda pair: (pair[0], -pair[1]))
            compat_seq = [row for row, _ in entries]
            wc_word = [letter for _, letter in entries]
            wc = cls.from_word_compatible(wc_word, compat_seq, length=length)
            assert wc.is_valid
            assert wc.perm == principal.perm
            wc_set.add(wc)
        return wc_set

    @cached_property
    def forest_weight(self):
        from schubmult.utils.tuple_utils import pad_tuple

        return pad_tuple(self.forest_invariant.forest.code, len(self))

    @property
    def grove_weight(self):
        return self.grove_invariant.length_vector

    @property
    def grove_invariant(self) -> tuple[int, ...]:
        r"""The composition of the grove containing ``self`` -- the inverse of
        :meth:`grove_wcs`.

        Every WCGraph of a fixed permutation lies in exactly one grove.  The grove
        is generated by the set-valued labelings of a single indexed forest ``F``,
        and its *base* member is the distinguished graph whose ``forest_weight``
        equals its ``length_vector`` (the minimal, one-element-per-node labeling).
        The grove weight is that base's ``forest_weight`` -- i.e. ``F.code`` padded
        to ``len(self)``.

        For a base member (``forest_weight == length_vector``) the grove weight
        coincides with the forest weight; for set-valued members it recovers the
        composition of the base rather than the (larger, occurrence-inflated)
        forest weight of ``self`` itself.
        """
        length = len(self)
        key = self.perm
        cached = WCGraph._grove_invariant_cache.get(key)
        if cached is None:
            mapping = {}
            graphs = WCGraph.all_wc_graphs(self.perm)
            for base in graphs:
                if base.forest_weight != base.length_vector:
                    continue

                for member in WCGraph.grove_wcs(base.forest_weight, len(base), base):
                    mapping[member] = base
            WCGraph._grove_invariant_cache[key] = mapping
        else:
            mapping = cached
        return mapping[self.normalize()].resize(length)

    @classmethod
    def one_row(cls, a: int) -> WCGraph:
        return cls([tuple(range(a, 0, -1))])

    @property
    @cache
    def omega_invariant(self):
        from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion

        word = list(reversed(self.perm_word))

        def word_to_pair_labeled(word):
            counts = {}
            out = []
            for a in word:
                aa = int(a)
                counts[aa] = counts.get(aa, 0) + 1
                out.append(letterpair(aa, counts[aa]))
            return tuple(out)

        return omega_insertion(word_to_pair_labeled(word))

    @property
    def forest_invariant(self):
        return self.omega_invariant[0]

    def flipped_co_wc(self):
        spet = self.resize(len(self.perm) - 1)
        refspet = []
        for i in range(1, spet.rows):
            refspet.append([])
            for j in range(1, spet.rows):
                if spet.has_element(i, j):
                    refspet[-1] = refspet.toggle_ref_at(spet.rows - i - j + 1, j)
        return refspet

    def to_reduced_compatible_set_sequence(self):
        word = []
        seq = self.compatible_sequence
        set_seq = []
        working_perm = Permutation([])
        for i, letter in enumerate(self.perm_word):
            if working_perm[letter - 1] > working_perm[letter]:
                for root_index in range(len(word)):
                    root = working_perm.right_root_at(root_index, word=word)
                    if root == (letter, letter + 1):
                        set_seq[root_index].add(seq[i])
                        break
            else:
                working_perm = working_perm.swap(letter - 1, letter)
                word.append(letter)
                set_seq.append({seq[i]})
        return tuple(word), tuple(tuple(sorted(s)) for s in set_seq)

    @classmethod
    def _from_root_dict(cls, root_dict: dict[tuple[int, int], set[int]], length: int | None = None) -> WCGraph:
        ret_word = []
        compat_seq = []
        while root_dict:
            max_pair = max([(root, v) for root, v in root_dict.items() if root[1] == root[0] + 1], key=lambda x: (max(x[1]), -x[0][1], x[0][0]))
            letter = max_pair[0][0]
            ret_word = [letter, *ret_word]

            compat_seq += [max(max_pair[1])]
            root_dict[max_pair[0]].remove(max(max_pair[1]))
            if len(root_dict[max_pair[0]]) == 0:
                root_dict = {Permutation([]).swap(letter - 1, letter).act_root(*root): v for root, v in root_dict.items() if root != max_pair[0]}

        return cls.from_word_compatible(ret_word, sorted(compat_seq), length=length)

    @classmethod
    def from_reduced_compatible_set_sequence(cls, word, set_seq, length=None):
        if len(word) != len(set_seq):
            raise ValueError("Word and set sequence must have the same length")
        working_set_seq = [set(s) for s in set_seq]
        root_dict = {}
        for i in range(len(word)):
            root_dict[Permutation._right_root_at(i, word=word)] = working_set_seq[i]

        return cls._from_root_dict(root_dict, length=length)

    def _snap_reduced(self):
        from .rc_graph import RCGraph

        red_word, seq = self.to_reduced_compatible_set_sequence()
        compat_seq = [min(v) for v in seq]
        return RCGraph.from_reduced_compatible(red_word, compat_seq, length=len(self))

    # def to_word_compatible_set_sequence(self):
    #     from .nilplactic import NilPlactic
    #     increasing, recording = self.hecke_invariant

    #     hw_recording, raise_seq = recording.to_highest_weight(length=len(self))

    #     hw_rc = NilPlactic.from_word(tuple(reversed(increasing.row_word))).hw_rc(len(self))
    #     word = hw_rc.perm_word
    #     hw_compat = []
    #     for i in range(len(hw_recording)):
    #         hw_compat.extend(hw_recording[i])
    #     hw_wc = WCGraph.from_word

    def _snap_min(self):
        from .increasing_tableau import IncreasingTableau
        from .set_valued_tableau import SetValuedTableau
        increasing_tab, recording_tab = self.hecke_invariant
        recording_tab = SetValuedTableau({coord: (min(labels),) for coord, labels in recording_tab.cells.items()})
        return WCGraph.from_word_compatible(*tuple(reversed(IncreasingTableau.hecke_column_uninsert_rsk(increasing_tab, recording_tab))))

    @classmethod
    def from_word_compatible(cls, word, seq, length=None):
        if length is None:
            length = max(seq, default=0)
        if length < max(seq, default=0):
            raise ValueError("Length must be at least the maximum value in the sequence")
        if len(word) != len(seq):
            raise ValueError("Word and sequence must have the same length")
        if any(seq[i] > seq[i + 1] for i in range(len(seq) - 1)):
            raise ValueError("Sequence must be weakly increasing")
        if any(word[i] <= word[i + 1] and seq[i] == seq[i + 1] for i in range(len(seq) - 1)):
            raise ValueError(f"{seq} is not compatible with {word}")
        result = [[] for _ in range(length)]
        for i in range(len(word)):
            result[seq[i] - 1] = [*result[seq[i] - 1], word[i]]
        return cls(tuple([tuple(row) for row in result]))

    def left_to_right_inversion_coords(self, index: int) -> tuple[int, int]:
        if index < 0 or index >= len(self.perm_word):
            raise ValueError(f"Index {index} out of range {self.perm.inv}")
        index_find = 0
        for i in range(len(self)):
            index_find2 = index_find + len(self[i])
            if index_find2 > index:
                break
            index_find = index_find2

        return (i + 1, self[i][(index - index_find)] - i)

    def left_to_right_inversion(self, index: int) -> tuple[int, int]:
        coords = self.left_to_right_inversion_coords(index)
        return self.right_root_at(*coords)

    def left_to_right_hecke_inversion(self, index: int) -> tuple[int, int]:
        coords = self.left_to_right_inversion_coords(index)
        return self.right_hecke_root_at(*coords)

    def right_root_at(self, i: int, j: int) -> tuple[int, int]:
        if i <= 0 or j <= 0:
            raise IndexError("i and j must be positive")
        if len(self.perm_word) > 0:
            index = self.bisect_left_coords_index(i, j)
            if index < len(self.perm_word):
                if self.left_to_right_inversion_coords(index) == (i, j):
                    return self.perm.right_root_at(index, word=self.perm_word)
                word_piece = list(self.perm_word[index:])
            else:
                word_piece = []
            refl = ~Permutation.ref_product(*word_piece)
            result = refl.act_root(i + j - 1, i + j)

        else:
            result = (i + j - 1, i + j)
        return result

    def right_hecke_root_at(self, i: int, j: int) -> tuple[int, int]:
        if i <= 0 or j <= 0:
            raise IndexError("i and j must be positive")
        if len(self.perm_word) > 0:
            index = self.bisect_left_coords_index(i, j)
            if index < len(self.perm_word):
                if self.left_to_right_inversion_coords(index) == (i, j):
                    return self.perm.right_hecke_root_at(index, word=self.perm_word)
                word_piece = list(self.perm_word[index:])
            else:
                word_piece = []
            refl = ~Permutation.hecke_ref_product(*word_piece)
            result = refl.act_root(i + j - 1, i + j)

        else:
            result = (i + j - 1, i + j)
        return result

    def polyvalue(self, x: Sequence[Expr], y: Sequence[Expr] | None = None, *, beta: Expr = None, prop_beta: bool = False, crystal: bool = False) -> Expr:
        if crystal:
            raise NotImplementedError("WCGraph is not a crystal graph")
        ret = S.One
        if beta is None:
            for i, row in enumerate(self):
                if y is None:
                    ret *= x[i + 1] ** len(row)
                else:
                    ret *= prod([x[i + 1] - y[row[j] - i] for j in range(len(row))])
            return ret
        for i, row in enumerate(self):
            if y is None:
                ret *= x[i + 1] ** len(row)
            else:
                ret *= prod([x[i + 1] + y[row[j] - i] + beta * x[i + 1] * y[row[j] - i] for j in range(len(row))])
        if prop_beta:
            ret *= beta ** (len(self.perm_word) - self.hecke_perm.inv)
        else:
            ret *= beta ** (len(self.perm_word))
        return ret

    @property
    def is_elem_sym(self) -> bool:
        cd = self.perm.trimcode
        if len(cd) == 0:
            return True
        if len(cd) == 1 and cd[0] == 1:
            return True
        spot = max([i for i, c in enumerate(cd) if c == 0], default=-1) + 1
        return all(c == 0 for c in cd[:spot]) and all(c == 1 for c in cd[spot:])
        #return False

    def crystal_length(self) -> int:
        return len(self)

    @property
    def hecke_invariant(self):
        from .increasing_tableau import IncreasingTableau
        return IncreasingTableau.hecke_column_insert_rsk(self.compatible_sequence, self.perm_word)

    @property
    def strong_hecke_invariant(self):
        """K-theoretic rectification of the diagonal-strip increasing tableau.

        Lay the ``perm_word`` out along a single anti-diagonal strip, reading
        from bottom-left to top-right, as a skew
        :class:`~schubmult.combinatorics.increasing_tableau.IncreasingTableau`,
        then repeatedly apply the simultaneous inner-corner K-theoretic down
        slide until the tableau is rectified (no inner corners remain). The
        rectified increasing tableau is the strong Hecke invariant.
        """
        from .increasing_tableau import IncreasingTableau

        word = self.perm_word
        m = len(word)
        if m == 0:
            return IncreasingTableau()
        # Cell (m-1-i, i) carries word[i]: word[0] sits bottom-left, word[-1]
        # top-right, so each letter occupies its own anti-diagonal cell.
        rows = tuple(tuple([0] * (m - 1 - r) + [word[m - 1 - r]]) for r in range(m))
        inner_shape = tuple(m - 1 - r for r in range(m))
        tab = IncreasingTableau(rows, inner_shape=inner_shape)
        while list(tab.iter_inner_corners):
            tab = tab.down_jdt_slide_all_inner_corners()
        return tab

    @cache
    def _convert_elem_rc(self):
        if self.is_elem_sym:
            from .rc_graph import RCGraph
            return next(iter(RCGraph.elem_sym_rcs(len(self.perm_word), self.perm.max_descent, weight=self.length_vector)))
        raise ValueError("Cannot convert non-elementary symmetric WCGraph to RCGraph")

    def elem_sym_wcs(p, k, weight=None):
        from schubmult import uncode
        perm = uncode([0] * (k - p) + [1] * p)
        return WCGraph.all_wc_graphs(perm, k, weight=weight)

    def to_rc_pieri(self):
        """Map this WCGraph to an RCGraph via ``_snap_reduced`` + Pieri insertion.

        Snap to the underlying reduced RCGraph, then reinsert the missing
        (excess) letters row-by-row: the number of letters missing from row
        ``i + 1`` is ``self.length_vector[i] - reduced.length_vector[i]``, and
        those rows (with multiplicity) are Pieri-inserted at ``perm.max_descent``.
        """
        from .rc_graph import RCGraph
        result = RCGraph([()]).resize(len(self))
        for i in range(len(self)):
            for j in range(len(self[i]) - 1, -1, -1):
                result = result.pieri_insert(self[i][j], [i + 1])
        return result

    @property
    def crystal_weight(self) -> tuple[int, ...]:
        return self.length_vector

    @property
    def excess(self):
        return len(self.perm_word) - self.perm.inv


    def raising_operator(self, i: int) -> WCGraph | None:
        from .increasing_tableau import IncreasingTableau

        if i <= 0 or i >= len(self):
            return None
        p_tab, q_tab = self.hecke_invariant
        q_next = q_tab.raising_operator(i)
        if q_next is None:
            return None
        if self._double:
            q_next = q_next.raising_operator(i)
            if q_next is None:
                return None
        try:
            compat, word = IncreasingTableau.hecke_column_uninsert_rsk(p_tab, q_next)
            result = WCGraph.from_word_compatible(word, compat, length=len(self))
        except Exception:
            return None
        if result.perm != self.perm or not result.is_valid:
            return None
        return result

    def lowering_operator(self, i: int) -> WCGraph | None:
        from .increasing_tableau import IncreasingTableau

        if i <= 0 or i >= len(self):
            return None
        p_tab, q_tab = self.hecke_invariant
        q_next = q_tab.lowering_operator(i)
        if q_next is None:
            return None
        if self._double:
            q_next = q_next.lowering_operator(i)
            if q_next is None:
                return None
        try:
            compat, word = IncreasingTableau.hecke_column_uninsert_rsk(p_tab, q_next)
            result = WCGraph.from_word_compatible(word, compat, length=len(self))
        except Exception:
            return None
        if result.perm != self.perm or not result.is_valid:
            return None
        return result

    @classmethod
    @cache
    def _extremal_weight(cls, perm: Permutation, length, h_inv):
        return min(wc.length_vector for wc in WCGraph.all_wc_graphs(perm, length) if wc.strong_hecke_invariant == h_inv)

    @cached_property
    def sorted_length_vector(self):
        return tuple(sorted(self.length_vector, reverse=True))


    @property
    def extremal_weight(self):
        # from schubmult.utils.tuple_utils import pad_tuple
        from schubmult.utils.tuple_utils import pad_tuple
        return pad_tuple(WCGraph._extremal_weight(self.perm, len(self), self.strong_hecke_invariant), len(self))

    @classmethod
    @cache
    def groth_to_schub(cls, groth_perm: Permutation, beta: Expr):
        from .pipe_dream import PipeDream

        boip = WCGraph.all_wc_graphs(groth_perm, len(groth_perm))
        ret = {}
        for wc in boip:
            cpdb = PipeDream.from_rc_graph(wc).co_pipe_dream()
            if cpdb.is_reduced:
                w0 = Permutation.w0(cpdb.rows)
                perm = cpdb.perm * w0
                ret[perm] = ret.get(perm, 0) + beta ** (perm.inv - groth_perm.inv)
        return ret

    @classmethod
    @cache
    def schub_to_groth(cls, schub_perm: Permutation, beta: Expr):
        from .pipe_dream import PipeDream
        from .rc_graph import RCGraph

        boip = RCGraph.all_rc_graphs(schub_perm, len(schub_perm.trimcode))
        ret = {}
        for rc in boip:
            cpdb = PipeDream.from_rc_graph(rc).co_pipe_dream()
            w0 = Permutation.w0(cpdb.rows)
            perm = cpdb.perm * w0
            ret[perm] = ret.get(perm, 0) + (-beta) ** (perm.inv - schub_perm.inv)
        return ret

    @classmethod
    @cache
    def grove_to_forest(cls, comp, beta: Expr):
        from .pipe_dream import PipeDream

        boip = WCGraph.grove_wcs(comp, len(comp))
        ret = {}
        for wc in boip:
            cpdb = PipeDream.from_rc_graph(wc).co_pipe_dream()
            if cpdb.is_reduced:
                rcc = cpdb.to_rc_graph()
                fweight = tuple((rcc.perm * Permutation.w0(wc.rows)).pad_code(len(comp)))
                ret[fweight] = ret.get(fweight, 0) + beta ** (sum(fweight) - sum(comp))
        return ret

    @cache
    def bisect_left_coords_index(self, row: int, col: int, lo: int = 0, hi: int | None = None) -> int:
        from bisect import bisect_left, bisect_right  # noqa: F401

        if hi is None:
            hi = len(self.perm_word)

        while lo < hi:
            mid = (lo + hi) // 2
            i, j = self.left_to_right_inversion_coords(mid)
            if i < row or (i == row and j > col):
                lo = mid + 1
            else:
                hi = mid
        return lo

    def vertical_cut(self, row: int) -> tuple[WCGraph, WCGraph]:
        if row < 0:
            raise ValueError("Row out of range")
        if row == 0:
            return self._rebuild(()), self
        if row == len(self):
            return self, self._rebuild(())
        if row >= len(self):
            raise ValueError("Row out of range")
        front = self._rebuild([*self[:row]])
        front = front.extend(max(len(self), len(front.perm.trimcode)) - row)
        flen = len(front)
        for _ in range(flen - row):
            front = front.zero_out_last_row()
        if row == len(self):
            back = self._rebuild([])
        else:
            back = self.rowrange(row, len(self))
        return (front, back)

    def disjoint_union(self, rc: WCGraph) -> WCGraph:
        if len(self) != len(rc):
            raise ValueError(f"{type(self).__name__}s must have the same number of rows")
        if self.perm.inv == 0:
            return rc
        rowmax = [max(self[i], default=0) for i in range(len(self))]
        N = max(rowmax) + 1
        shift_rc = self._rebuild([tuple([a + N for a in row]) for row in rc]).resize(len(rc) + N)
        rc_self = self.resize(len(rc) + N)
        return self._rebuild([shift_rc[i] + rc_self[i] for i in range(len(rc_self))])

    @property
    def is_quasi_yamanouchi(self) -> bool:
        for i in range(1, len(self)):
            if max(self[i], default=0) <= min(self[i - 1], default=0) and min(self[i - 1], default=0) >= i + 1:
                return False
        return True

    @property
    def dst(self):
        if self.is_quasi_yamanouchi:
            return self
        rows = [*self]
        changed = True
        while changed:
            changed = False
            for i in range(1, len(rows)):
                if max(rows[i], default=0) <= min(rows[i - 1], default=0) and min(rows[i - 1], default=0) >= i + 1:
                    rows[i] = sorted({*rows[i - 1], *rows[i]}, reverse=True)
                    rows[i - 1] = []
                    changed = True
                    break
        ret = WCGraph([tuple(row) for row in rows])
        if not ret.is_quasi_yamanouchi:
            raise ValueError("DST failed to produce a quasi-Yamanouchi graph")
        return ret

    @cache
    def squash_product(self, rc: WCGraph) -> WCGraph:
        combined_rc = self.disjoint_union(rc)
        while len(combined_rc) > len(self):
            combined_rc = combined_rc.zero_out_last_row()
        return combined_rc

    @cache
    def zero_out_last_row(self) -> WCGraph:
        r"""Zero out the (empty) last row, realizing the Weigandt/Lascoux
        transition (``writing/weigandt_bumpless.tex``, Theorem "transition").

        The graph is transported to a marked bumpless pipedream, the underlying
        BPD is resized down by one row (dropping the maximal-corner row), and the
        result is transported back.  Weight is preserved; the permutation changes
        according to the transition.  Unlike the previous Hecke-insertion route,
        this is total: it works for non-core-reduced graphs as well.

        Cached: the full MBPD round trip is memoized on the (hashable) graph,
        so repeated ``squash_product`` calls reuse the result."""
        if len(self) == 0:
            return self
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")
        rows = len(self) - 1
        reduced = self.to_mbpd().zero_out_last_row(rows)
        return WCGraph.from_mbpd(reduced).resize(rows)

    _z_cache: ClassVar[dict[WCGraph, set[WCGraph]]] = {}

    def right_zero_act(self) -> set[WCGraph]:
        from schubmult import AGx, uncode
        if self.perm.inv == 0:
            return {self._rebuild([*self, ()])}

        if self in WCGraph._z_cache:
            return WCGraph._z_cache[self]

        up_perms = {perm for (perm, _), v in (AGx(self.perm, len(self)) * AGx(uncode([0]), 1)).items() if v != 0}
        # stack =
        # ci inverse times descent

        # def _cycle_comb(pivot_set):
        #     pivot_list = sorted(pivot_set, reverse=True)
        #     #maxd, b = self.maximal_corner
        #     maxd = self.perm.max_descent + 1
        #     cycle = [maxd, *pivot_list]
        #     cycle_arr = list(range(1, len(self.perm) + 1))
        #     for i in range(len(cycle) - 1):
        #         cycle_arr[cycle[i] - 1] = cycle[i + 1]
        #     cycle_arr[cycle[-1] - 1] = cycle[0]
        #     cycle_perm = Permutation(cycle_arr)
        #     return cycle_perm

        # for r in range(1, self.perm.max_descent):
        #     for comb in itertools.combinations(range(1, self.perm.max_descent), r):
        #         cyc = _cycle_comb(comb)
        #         newperm = self.perm * (~cyc)
        #         if newperm.inv == self.perm.inv - r:
        #             a = self.perm.max_descent + 1
        #             for b_prime in range(a + 1, len(self.perm) + 1):
        #                 if newperm[a - 1] < newperm[b_prime - 1]:
        #                     new_new_perm = newperm.swap(a - 1, b_prime - 1)
        #                     if new_new_perm.max_descent == self.perm.max_descent + 1 and new_new_perm.inv == self.perm.inv - r + 1:
        #                         up_perms.add(new_new_perm)
        #             #up_perms.add(newperm)
        # #{perm for perm in Permutation.all_permutations(len(self.perm) + 1) if self.perm == perm for self.perm in pivot_transition(perm)}

        rc_set = set()

        for perm in up_perms:
            for wc in type(self).all_wc_graphs(perm, len(self) + 1, weight=(*self.length_vector, 0)):
                if wc.zero_out_last_row() == self:
                    rc_set.add(wc)

        WCGraph._z_cache[self] = rc_set
        return rc_set

    @classmethod
    def principal_wc(cls, perm, length):
        from schubmult.combinatorics.rc_graph import RCGraph
        return cls(RCGraph.principal_rc(perm, length))

    def product(self, other: SchubertMonomialGraph) -> dict[WCGraph, int]:
        """Compute the product of this WC graph with another."""
        from schubmult.utils.perm_utils import add_perm_dict
        self_len = len(self)
        if self.perm.inv == 0:
            return {self._rebuild([*self, *other.shiftup(self_len)]): 1}
        num_zeros = max(len(other), len(other.perm))
        assert len(self.perm.trimcode) <= self_len, f"{self=}, {self.perm=}"
        base_rc = self
        buildup_module = {base_rc: 1}

        for _ in range(num_zeros):
            new_buildup_module = {}
            for rc, coeff in buildup_module.items():
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(rc.right_zero_act(), coeff))
            buildup_module = new_buildup_module
        ret_module = {}
        other_shifted = other.shiftup(self_len)
        target_len = self_len + len(other)

        for rc, coeff in buildup_module.items():
            new_rc = type(rc)([*rc[:self_len], *other_shifted])
            assert len(new_rc) == target_len
            if new_rc.is_valid and len(new_rc.perm.trimcode) <= len(new_rc):
                ret_module = add_perm_dict(ret_module, {new_rc: coeff})

        return ret_module

    def _upieri_insert_row(self, row, descent, dict_by_a, dict_by_b, num_times, start_index=-1, backwards=True, reflection_rows=None, target_row=None, left=False):
        working_rc = self
        if descent is not None and row > descent:
            raise ValueError("All rows must be less than or equal to descent")

        i = start_index
        new_reflections = []  # Track reflections added in THIS call

        if i == -1:
            if backwards or (not backwards and left):
                i = 0
            else:
                i = max(self[row - 1], default=0) + descent + 5
        num_done = 0
        flag = True
        attempts_without_progress = 0
        last_num_done = -1
        max_attempts = 100
        while num_done < num_times:
            if num_done == last_num_done:
                attempts_without_progress += 1
            last_num_done = num_done
            if attempts_without_progress > max_attempts:
                raise ValueError(
                    f"Pieri insertion made no progress after {max_attempts} attempts ({row=}, {descent=}, {num_times=}, {left=}, {backwards=})",
                )
            if i <= 1 and (not backwards or (backwards and left)):
                i = working_rc.cols + descent + 5
            if not backwards or (backwards and left):
                i -= 1
            else:
                i += 1
            flag = False

            if not working_rc.has_element(row, i):
                if left:
                    a, b = working_rc.left_root_at(row, i)
                else:
                    a, b = working_rc.right_root_at(row, i)
                if a < b:
                    flag = False
                    if _is_row_root(descent, (a, b)) and b not in dict_by_b:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[a] = dict_by_a.get(a, set())
                        dict_by_a[a].add(b)
                        dict_by_b[b] = a
                        if reflection_rows is not None and target_row is not None:
                            reflection_rows[(a, b)] = target_row
                            new_reflections.append((a, b))
                        flag = True
                    elif a in dict_by_b and b > descent and b not in dict_by_b:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        dict_by_a[dict_by_b[a]].add(b)
                        dict_by_b[b] = dict_by_b[a]
                        if reflection_rows is not None and target_row is not None:
                            reflection_rows[(dict_by_b[a], b)] = target_row
                            new_reflections.append((dict_by_b[a], b))
                        flag = True
                    elif descent is None:
                        working_rc = working_rc.toggle_ref_at(row, i)
                        flag = True
                if flag:
                    num_done += 1
                    attempts_without_progress = 0
                if not left and row > 1 and not working_rc.is_valid:
                    working_rc = working_rc._upieri_rectify(row - 1, descent, dict_by_a, dict_by_b, backwards=backwards, reflection_rows=reflection_rows, left=left)  # minus one?
                if left and row < len(working_rc) and not working_rc.is_valid:
                    working_rc = working_rc._upieri_rectify(row + 1, descent, dict_by_a, dict_by_b, backwards=backwards, reflection_rows=reflection_rows, left=left)
        return working_rc, new_reflections

    def _upieri_rectify(self, row_below, descent, dict_by_a, dict_by_b, backwards=True, reflection_rows=None, target_row=None, left=False):
        working_rc = self
        if working_rc.is_valid:
            return working_rc
        if row_below == 0 and not left:
            assert working_rc.is_valid, f"{working_rc=}, {dict_by_a=}, {dict_by_b=}"
            return working_rc
        if left and row_below > len(working_rc):
            if working_rc.is_valid:
                return working_rc
            raise ValueError(f"Left upieri rectify exhausted rows without reaching validity ({row_below=}, {len(working_rc)=})")
        extra = descent if descent is not None else 0
        the_range = range(1, working_rc.cols + extra + 5)
        if left:
            the_range = reversed(the_range)
        for j in the_range:
            flag = False
            if working_rc.is_valid:
                return working_rc

            if working_rc.has_element(row_below, j):
                if left:
                    a, b = working_rc.left_root_at(row_below, j)
                else:
                    a, b = working_rc.right_root_at(row_below, j)

                top, bottom = max(a, b), min(a, b)

                if a < b:
                    continue

                if bottom in dict_by_a and top in dict_by_a[bottom]:
                    new_rc = working_rc.toggle_ref_at(row_below, j)
                    dict_by_a[bottom].remove(top)
                    if len(dict_by_a[bottom]) == 0:
                        del dict_by_a[bottom]
                    del dict_by_b[top]
                    working_rc = new_rc
                    flag = True

                elif bottom in dict_by_b and top in dict_by_b and dict_by_b[top] == dict_by_b[bottom]:
                    new_rc = working_rc.toggle_ref_at(row_below, j)
                    dict_by_a[dict_by_b[bottom]].remove(top)

                    if len(dict_by_a[dict_by_b[top]]) == 0:
                        del dict_by_a[dict_by_b[top]]
                    del dict_by_b[top]
                    flag = True
                    working_rc = new_rc
                elif descent is None:
                    working_rc = working_rc.toggle_ref_at(row_below, j)
                    flag = True
                else:
                    raise ValueError(f"Could not rectify at {(row_below, j)} with root {(a, b)}")
                if flag:
                    working_rc, _ = working_rc._upieri_insert_row(
                        row_below,
                        descent,
                        dict_by_a,
                        dict_by_b,
                        num_times=1,
                        backwards=backwards,
                        reflection_rows=reflection_rows,
                        target_row=target_row,
                        left=left,
                    )
        if left:
            return working_rc._upieri_rectify(row_below + 1, descent, dict_by_a, dict_by_b, backwards=backwards, reflection_rows=reflection_rows, target_row=target_row, left=left)
        return working_rc._upieri_rectify(row_below - 1, descent, dict_by_a, dict_by_b, backwards=backwards, reflection_rows=reflection_rows, target_row=target_row, left=left)

    # VERIFY
    def upieri_insert(self, descent, rows, return_reflections=False, backwards=True, left=False):
        dict_by_a = {}
        dict_by_b = {}
        reflection_rows = {}  # Track which row each reflection was added to
        # row is descent
        # inserting times

        working_rc = self._rebuild([*self])
        if len(rows) == 0:
            if return_reflections:
                return working_rc, ()
            return self
        rows_grouping = {}

        for r in rows:
            rows_grouping[r] = rows_grouping.get(r, 0) + 1
        if max(rows) > len(working_rc):
            working_rc = working_rc.extend(max(rows) - len(working_rc))
        rows = sorted(rows, reverse=not left)
        reflections = []
        for row in sorted(rows_grouping.keys(), reverse=not left):
            num_times = rows_grouping[row]
            last_working_rc = working_rc
            working_rc, new_reflections = working_rc._upieri_insert_row(row, descent, dict_by_a, dict_by_b, num_times, backwards=backwards, reflection_rows=reflection_rows, target_row=row, left=left)
            reflections += new_reflections
            if (not left and row > 1) and not working_rc.is_valid:
                working_rc = working_rc._upieri_rectify(row - 1, descent, dict_by_a, dict_by_b, backwards=backwards, reflection_rows=reflection_rows, target_row=row, left=left)  # minus one?
            elif (left and row < len(working_rc)) and not working_rc.is_valid:
                working_rc = working_rc._upieri_rectify(row + 1, descent, dict_by_a, dict_by_b, backwards=backwards, reflection_rows=reflection_rows, target_row=row, left=left)
            try:
                assert len(working_rc[row - 1]) == len(last_working_rc[row - 1]) + num_times
            except AssertionError:
                raise
        if return_reflections:
            # Build list of (row, reflection) pairs
            return working_rc, tuple(reflections)
        return working_rc

    @staticmethod
    def _strict_decreasing_reflection_rows(max_reflection: int, k: int) -> tuple[tuple[int, ...], ...]:
        if k < 0:
            raise ValueError("k must be nonnegative")
        if k == 0:
            return ((),)
        if max_reflection < k:
            return ()
        # Choose k reflections and store each row in strict decreasing order,
        # matching the RC/WC internal row convention.
        return tuple(tuple(sorted(combo, reverse=True)) for combo in itertools.combinations(range(1, max_reflection + 1), k))

    @staticmethod
    @cache
    def _solve_shifted_right_hecke_factors(w_prime: Permutation, w: Permutation) -> tuple[Permutation, ...]:
        n = len(w)
        if n == 0:
            return (Permutation([]),)
        factors = []
        for wpp in Permutation.all_permutations(n - 1):
            if w_prime @ wpp.shiftup(1) == w:
                factors.append(wpp)
        return tuple(factors)

    @classmethod
    def pull_out_var_hecke(cls, w: Permutation, k: int) -> set[tuple[tuple[int, ...], Permutation]]:
        """Hecke/Demazure analogue of pull_out_var(1, w) for fixed first-row size.

        Returns all pairs ``(row, wpp)`` with:
        - ``row`` a strictly decreasing tuple of simple reflections of size ``k``
        - ``wpp`` a permutation such that
          ``Permutation.ref_product(*row) @ wpp.shiftup(1) == w``.
        """
        w = Permutation(w)
        if k < 0:
            raise ValueError("k must be nonnegative")
        if len(w) == 0:
            if k == 0:
                return {((), Permutation([]))}
            return set()

        ret: set[tuple[tuple[int, ...], Permutation]] = set()
        for row in cls._strict_decreasing_reflection_rows(len(w) - 1, k):
            w_prime = Permutation.ref_product(*row)
            for wpp in cls._solve_shifted_right_hecke_factors(w_prime, w):
                ret.add((row, wpp))
        return ret

    _graph_cache: dict[tuple[Permutation, int], set[WCGraph]] = {}  # noqa: RUF012
    _cache_by_weight: dict[tuple[Permutation, tuple[int, ...]], set[WCGraph]] = {}  # noqa: RUF012
    _grove_invariant_cache: dict[tuple[Permutation, int], dict[WCGraph, tuple[int, ...]]] = {}  # noqa: RUF012

    _slow_graph_cache: dict[tuple[Permutation, int], set[WCGraph]] = {}  # noqa: RUF012
    _slow_cache_by_weight: dict[tuple[Permutation, tuple[int, ...]], set[WCGraph]] = {}  # noqa: RUF012

    @classmethod
    def all_wc_graphs_slow(cls, perm: Permutation, length: int = -1, weight: tuple[int, ...] | None = None) -> set[WCGraph]:
        """Generate all WC graphs with Demazure permutation ``perm`` and fixed row count.

        The recursion mirrors ``all_rc_graphs`` but uses Hecke/Demazure pull-out on the
        first row via ``pull_out_var_hecke``.
        """
        perm = Permutation(perm)
        if length < 0:
            # Non-empty rows can only occur in indices 1..len(perm)-1.
            length = perm.max_descent
        if length < perm.max_descent:
            return set()
        if weight is not None and len(weight) != length:
            raise ValueError("Weight must have length equal to the number of rows")

        if weight is not None:
            if (perm, tuple(weight)) in cls._slow_cache_by_weight:
                return cls._slow_cache_by_weight[(perm, tuple(weight))]
        elif (perm, length) in cls._slow_graph_cache:
            return cls._slow_graph_cache[(perm, length)]

        if perm.inv == 0:
            if weight is not None and sum(weight) != 0:
                cls._slow_cache_by_weight[(perm, tuple(weight))] = set()
                return set()
            ret = {cls([()] * length)}
            if weight is not None:
                cls._slow_cache_by_weight[(perm, tuple(weight))] = ret
            else:
                cls._slow_graph_cache[(perm, length)] = ret
            return ret

        ret: set[WCGraph] = set()
        first_row_sizes = [weight[0]] if weight is not None else list(range(len(perm)))
        for k in first_row_sizes:
            for row, reduced_perm in cls.pull_out_var_hecke(perm, k):
                old_set = cls.all_wc_graphs_slow(reduced_perm, length=length - 1, weight=(None if weight is None else weight[1:]))
                for old_wc in old_set:
                    new_rows = [row, *[tuple(a + 1 for a in rr) for rr in old_wc]]
                    nwc = cls(new_rows)
                    if Permutation.ref_product(*row) @ reduced_perm.shiftup(1) != perm:
                        continue
                    if len(nwc) != length:
                        continue
                    if weight is not None and nwc.length_vector[0] != weight[0]:
                        continue
                    ret.add(nwc)

        if weight is not None:
            cls._slow_cache_by_weight[(perm, tuple(weight))] = ret
        else:
            cls._slow_graph_cache[(perm, length)] = ret
        return ret

    @classmethod
    def all_wc_graphs(cls, perm: Permutation, length: int | None = None, weight: tuple[int, ...] | None = None, *, check_length=False, do_cache=True) -> set[WCGraph]:
        from schubmult.utils.schub_lib import pull_out_var
        # print(f"{weight=}")
        if check_length and length is not None and (length > 0 and length < len(perm.trimcode)):
            raise ValueError(f"Length must be at least the last descent of the permutation, permutation has {len(perm.trimcode)} rows and {perm=}, got {length=}")
        if length is None:
            length = len(perm.trimcode)
        if length < 0:
            raise ValueError("Length must be nonnegative")
        if length == 0:
            if perm.inv == 0:
                return {cls([])}
            return set()

        if weight is not None and len(weight) != length:
            raise ValueError("Weight must have length equal to the number of rows")
        if weight is not None and do_cache:
            if (perm, tuple(weight)) in cls._cache_by_weight:
                return cls._cache_by_weight[(perm, tuple(weight))]
        elif do_cache and (perm, length) in cls._graph_cache:
            return cls._graph_cache[(perm, length)]
        # if perm.inv == 0:
        #     if do_cache:
        #         cls._graph_cache[(perm, length)] = {cls([()] * length if length > 0 else [])}
        #         return cls._graph_cache[(perm, length)]
        #     return {cls([()] * length if length > 0 else [])}
        ret = set()
        pm = perm
        L = [*pull_out_var(1, pm)]
        if length == 1:
            if (weight is not None and weight[0] != pm.inv) or pm.max_descent > 1:
                return set()
            ret = {WCGraph.one_row(pm.inv)}
        else:
            for _, new_perm in L:
                new_row_base = sorted([new_perm[i] for i in range(max(len(pm), len(new_perm))) if new_perm[i] == pm[i + 1]], reverse=True)
                new_perm_up = new_perm.shiftup(1)
                addable_descents = set()
                for dd in set(range(1, len(perm))) - set(new_row_base):
                    test_row = sorted((*new_row_base, dd), reverse=True)
                    if Permutation.ref_product(*test_row) @ new_perm_up == perm:
                        addable_descents.add(dd)
                stack = [(new_row_base, addable_descents)]
                if weight is not None:
                    oldset = cls.all_wc_graphs(new_perm, length=length - 1, weight=weight[1:], do_cache=do_cache)
                else:
                    oldset = cls.all_wc_graphs(new_perm, length=length - 1, do_cache=do_cache)
                all_new_rows = set()
                while len(stack) > 0:
                    new_row, addable_descents = stack.pop()
                    if weight is None or len(new_row) == weight[0]:
                        all_new_rows.add(tuple(new_row))
                    if weight is None or len(new_row) < weight[0]:
                        for d in addable_descents:
                            new_new_row = sorted((*new_row, d), reverse=True)
                            addable_descents2 = set()
                            for dd in set(range(1, len(perm))) - set(new_new_row):
                                test_row = sorted((*new_new_row, dd), reverse=True)
                                if Permutation.ref_product(*test_row) @ new_perm_up == perm:
                                    addable_descents2.add(dd)
                            stack.append((new_new_row, addable_descents2))
                for old_wc in oldset:
                    oldup = [tuple([row[i] + 1 for i in range(len(row))]) for row in old_wc]
                    for new_row in all_new_rows:
                        nwc = cls([tuple(new_row), *oldup])
                        assert nwc.perm == perm
                        assert len(nwc) == length, f"{nwc=}, {length=}, {old_wc=}, {new_row=}"
                        if weight is not None:
                            assert nwc.length_vector[1:] == tuple(weight[1:]), f"{nwc=}, {weight=}"
                        ret.add(nwc)
        if weight is not None and do_cache:
            cls._cache_by_weight[(perm, tuple(weight))] = ret
        else:
            if do_cache:
                cls._graph_cache[(perm, length)] = ret
        # print(f"{perm=} {weight=} {length=} {len(ret)=}")
        return ret

    @classmethod
    def grothendieck_polynomial_via_wc(cls, perm: Permutation, x: Sequence[Expr], beta: Expr, length: int = -1):
        """Compute a Grothendieck candidate by summing WC graph monomials.

        Uses weight monomials with ``beta`` exponent ``|word|-inv(perm)``.
        """
        ret = S.Zero
        for wc in cls.all_wc_graphs(Permutation(perm), length=length):
            ret += wc.polyvalue(x, beta=beta, prop_beta=True)
        return ret

    def __getitem__(self, key: int | tuple[int, int]) -> tuple[int, ...] | int:
        if isinstance(key, int):
            return tuple(self)[key]
        if isinstance(key, tuple):
            i, j = key
            if isinstance(i, slice):
                return tuple([self[a, j] for a in range(len(self))[i]])
            if isinstance(j, slice):
                return tuple([self[i, b] for b in range(self.cols)[j]])
            if not self.has_element(i + 1, j + 1):
                return None
            return i + j + 1
        is_slice = isinstance(key, slice)
        if is_slice:
            return tuple(tuple(self)[n] for n in range(len(self))[key])
        raise ValueError(f"Bad indexing {key=}")

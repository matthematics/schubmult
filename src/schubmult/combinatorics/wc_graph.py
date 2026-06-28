from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from functools import cache, cached_property
from itertools import combinations

from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.schubert_monomial_graph import SchubertMonomialGraph
from schubmult.symbolic import Expr, S, prod
from schubmult.utils._grid_print import GridPrint

#from schubmult.utils.perm_utils import _is_compatible

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
        if compat_seq[i-1] > compat_seq[i]:
            return False
        if word[i - 1] <= word[i] and compat_seq[i-1] == compat_seq[i]:
            return False
        if compat_seq[i] > word[i]:
            return False
    return True


class WCGraph(SchubertMonomialGraph, GridPrint, tuple):
    """Word-compatible graph.

    Internal representation matches RCGraph: a tuple of rows where row i (0-indexed)
    is a strictly decreasing tuple of positive integers >= i + 1.

    For WCGraph, the associated permutation is the Demazure product of the
    concatenated row word (1-indexed simple reflections).
    """

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
        return hash(tuple(self))

    def trans_co_pipe(self):
        #return self._rebuild([tuple(reversed([])) for i, row in enumerate(reversed(self))])
        new_wc = WCGraph([]).resize(2*len(self))
        for i in range(1,self.rows + 1):
            for j in range(1, self.cols + 1):
                if not self.has_element(i,j):
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

    @property
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
        return self.resize(len(self.perm.trimcode))

    def resize(self, new_length: int) -> WCGraph:
        if new_length < len(self):
            return self.rowrange(0, new_length)
        return self.extend(new_length - len(self))

    def rowrange(self, start: int, end: int | None = None) -> WCGraph:
        if not end:
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

    @property
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

    @property
    def compatible_sequence(self) -> tuple[int, ...]:
        seq = []
        for i in range(len(self)):
            for _ in range(len(self[i])):
                seq.append(i + 1)
        return tuple(seq)

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
            raise ValueError(f"{type(self).__name__}s must have at least as many rows")
        if self.perm.inv == 0:
            return rc
        rowmax = [max(self[i], default=0) for i in range(len(self))]
        N = max(rowmax)
        shift_rc = self._rebuild([tuple([a + N for a in row]) for row in rc]).resize(len(rc) + N)
        rc_self = self.resize(len(rc) + N)
        return self._rebuild([shift_rc[i] + rc_self[i] for i in range(len(rc_self))])

    def squash_product(self, rc: WCGraph) -> WCGraph:
        combined_rc = self.disjoint_union(rc)
        while len(combined_rc) > len(self):
            combined_rc = combined_rc.zero_out_last_row()
        return combined_rc

    @cache
    def zero_out_last_row(self) -> WCGraph:
        if len(self) == 0:
            return self
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")

        reduced = self
        if self.perm.inv != len(self.perm_word):
            reduction_working = self
            removed_positive_roots: Counter[tuple[int, int]] = Counter()
            changed = True
            while changed:
                changed = False
                for index in range(len(reduction_working.perm_word) - 1, -1, -1):
                    a, b = reduction_working.left_to_right_inversion(index)
                    if b < a:
                        matching_positive = (b, a)
                        delete_index = None
                        for index2 in range(index + 1, len(reduction_working.perm_word)):
                            if reduction_working.left_to_right_inversion(index2) == matching_positive:
                                delete_index = index2
                                break
                        if delete_index is None:
                            continue
                        removed_positive_roots[matching_positive] += 1
                        row, col = reduction_working.left_to_right_inversion_coords(delete_index)
                        reduction_working = reduction_working.toggle_ref_at(row, col)
                        changed = True
                        break

            reduced = self._rebuild(reduction_working._snap_reduced().zero_out_last_row())
        else:
            return self._rebuild(reduced._snap_reduced().zero_out_last_row())

        working = reduced
        root_map = {reduced.left_to_right_inversion(i): reduced.left_to_right_inversion(i) for i in range(reduced.perm.inv)}
        target_positive_roots: Counter[tuple[int, int]] = Counter()
        for root, mult in removed_positive_roots.items():
            mapped_root = root_map.get(root)
            if mapped_root is None:
                continue
            target_positive_roots[mapped_root] += mult

        target_row_weights = tuple(self.length_vector[:-1])
        if len(target_row_weights) != len(working):
            raise ValueError("Row count mismatch while rebuilding WCGraph after zero_out_last_row")

        changed = True
        while changed and sum(target_positive_roots.values()) > 0:
            changed = False
            max_col = max(working.cols + sum(target_positive_roots.values()) + 1, len(self.perm) + len(self.perm_word) + 8)
            for row_index in range(len(working)):
                if len(working[row_index]) >= target_row_weights[row_index]:
                    continue
                for col in range(1, max_col + 1):
                    if working.has_element(row_index + 1, col):
                        continue
                    candidate = working.toggle_ref_at(row_index + 1, col)

                    # Index of the inserted crossing in the candidate perm_word.
                    row_value = row_index + col
                    try:
                        pos_in_row = candidate[row_index].index(row_value)
                    except ValueError:
                        continue
                    ins_index = sum(len(candidate[r]) for r in range(row_index)) + pos_in_row

                    tail = candidate.perm_word[ins_index + 1 :]
                    base_root = (row_index + col, row_index + col + 1)
                    matched_root = None

                    # Allow temporary chopping of trailing reflections on the right.
                    for chop in range(len(tail) + 1):
                        piece = tail[: len(tail) - chop]
                        if len(piece) == 0:
                            root = base_root
                        else:
                            apply = ~Permutation.ref_product(*piece)
                            root = apply.act_root(*base_root)
                        if target_positive_roots[root] > 0:
                            matched_root = root
                            break

                    if matched_root is not None:
                        working = candidate
                        target_positive_roots[matched_root] -= 1
                        if target_positive_roots[matched_root] == 0:
                            del target_positive_roots[matched_root]
                        changed = True
                        break
                if changed:
                    break

        if working.length_vector != target_row_weights:
            red_word, exist_compat = working.to_reduced_compatible_set_sequence()
            diffs = []
            for i in range(len(working)):
                if len(working[i]) < target_row_weights[i]:
                    diffs.extend([i + 1] * (target_row_weights[i] - len(working[i])))
            working2 = working.resize(len(self))
            col = 1
            while working2.perm.max_descent != len(self):
                working2 = self._rebuild([*working, tuple(range(len(self) + col - 1, len(self) - 1, - 1))])
                col += 1
            working2 = working2.upieri_insert(len(self), diffs).resize(len(self) - 1).resize(len(self)).zero_out_last_row()
            return working2

        return working

    def right_zero_act(self) -> set[WCGraph]:
        raise NotImplementedError("right_zero_act is implemented in RCGraph")

    def product(self, other: SchubertMonomialGraph) -> dict[WCGraph, int]:
        raise NotImplementedError("product is implemented in RCGraph")

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
        return tuple(tuple(sorted(combo, reverse=True)) for combo in combinations(range(1, max_reflection + 1), k))

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

    @classmethod
    def all_wc_graphs(cls, perm: Permutation, length: int = -1, weight: tuple[int, ...] | None = None) -> set[WCGraph]:
        """Generate all WC graphs with Demazure permutation ``perm`` and fixed row count.

        The recursion mirrors ``all_rc_graphs`` but uses Hecke/Demazure pull-out on the
        first row via ``pull_out_var_hecke``.
        """
        perm = Permutation(perm)
        if length < 0:
            # Non-empty rows can only occur in indices 1..len(perm)-1.
            length = max(len(perm) - 1, 0)
        if weight is not None and len(weight) != length:
            raise ValueError("Weight must have length equal to the number of rows")

        if weight is not None:
            wkey = (perm, tuple(weight))
            if wkey in cls._cache_by_weight:
                return cls._cache_by_weight[wkey]
        else:
            key = (perm, length)
            if key in cls._graph_cache:
                return cls._graph_cache[key]

        if length == 0:
            ret = {cls(())} if perm.inv == 0 else set()
            if weight is not None:
                cls._cache_by_weight[(perm, tuple(weight))] = ret
            else:
                cls._graph_cache[(perm, length)] = ret
            return ret

        if perm.inv == 0:
            ret = {cls([()] * length)}
            if weight is not None:
                cls._cache_by_weight[(perm, tuple(weight))] = ret
            else:
                cls._graph_cache[(perm, length)] = ret
            return ret

        ret: set[WCGraph] = set()
        first_row_sizes = [weight[0]] if weight is not None else list(range(len(perm)))
        for k in first_row_sizes:
            for row, reduced_perm in cls.pull_out_var_hecke(perm, k):
                old_set = cls.all_wc_graphs(reduced_perm, length=length - 1, weight=(None if weight is None else weight[1:]))
                for old_wc in old_set:
                    new_rows = [row, *[tuple(a + 1 for a in rr) for rr in old_wc]]
                    nwc = cls(new_rows)
                    if nwc.perm != perm:
                        continue
                    if len(nwc) != length:
                        continue
                    if weight is not None and nwc.length_vector != tuple(weight):
                        continue
                    ret.add(nwc)

        if weight is not None:
            cls._cache_by_weight[(perm, tuple(weight))] = ret
        else:
            cls._graph_cache[(perm, length)] = ret
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

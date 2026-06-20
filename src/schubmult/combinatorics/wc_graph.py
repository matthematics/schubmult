from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from functools import cache, cached_property
from itertools import combinations

from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.schubert_monomial_graph import SchubertMonomialGraph
from schubmult.symbolic import Expr, S, prod
from schubmult.utils._grid_print import GridPrint


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
    def from_reduced_compatible_set_sequence(cls, word, set_seq, length=None):
        ret_word = []
        compat_seq = []
        if len(word) != len(set_seq):
            raise ValueError("Word and set sequence must have the same length")
        working_set_seq = [set(s) for s in set_seq]
        root_dict = {}
        for i in range(len(word)):
            root_dict[Permutation._right_root_at(i, word=word)] = working_set_seq[i]

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
        index = 0
        for row_index in range(i - 1):
            index += len(self[row_index])
        index += self[i - 1].index(i + j - 1)
        word_piece = list(self.perm_word[index + 1 :])
        if len(word_piece) == 0:
            return (i + j - 1, i + j)
        apply = ~Permutation.ref_product(*word_piece)
        ret = apply.act_root(i + j - 1, i + j)
        if ret[0] > ret[1]:
            ret = (ret[1], ret[0])
        return ret

    def right_hecke_root_at(self, i: int, j: int) -> tuple[int, int]:
        if i <= 0 or j <= 0:
            raise IndexError("i and j must be positive")
        index = 0
        for row_index in range(i - 1):
            index += len(self[row_index])
        index += self[i - 1].index(i + j - 1)
        word_piece = list(self.perm_word[index + 1 :])
        if len(word_piece) == 0:
            return (i + j - 1, i + j)
        apply = ~Permutation.hecke_ref_product(*word_piece)
        ret = apply.act_root(i + j - 1, i + j)
        if ret[0] > ret[1]:
            ret = (ret[1], ret[0])
        return ret

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
    def zero_out_last_row(self) -> WCGraph:
        if len(self) == 0:
            return self
        if len(self[-1]) != 0:
            raise ValueError("Last row not empty")

        from .rc_graph import RCGraph

        reduction_working = RCGraph(self)
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

        rc = reduction_working
        reduced = rc.zero_out_last_row()
        if self.perm.inv == len(self.perm_word):
            return WCGraph(reduced)

        working = WCGraph(reduced)
        root_map = {rc.left_to_right_inversion(i): reduced.left_to_right_inversion(i) for i in range(rc.perm.inv)}
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

        if tuple(len(row) for row in working) != target_row_weights:
            raise ValueError(
                f"Failed to preserve row weights in WCGraph.zero_out_last_row: expected {target_row_weights}, got {tuple(len(row) for row in working)}",
            )
        if len(working.perm_word) != len(self.perm_word):
            raise ValueError(
                f"Failed to preserve crossing count in WCGraph.zero_out_last_row: expected {len(self.perm_word)}, got {len(working.perm_word)}",
            )
        return working

    def right_zero_act(self) -> set[WCGraph]:
        raise NotImplementedError("right_zero_act is implemented in RCGraph")

    def product(self, other: SchubertMonomialGraph) -> dict[WCGraph, int]:
        raise NotImplementedError("product is implemented in RCGraph")

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

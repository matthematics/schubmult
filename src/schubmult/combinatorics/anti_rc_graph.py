from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from functools import cached_property

import numpy as np

from schubmult.combinatorics.crystal_graph import CrystalGraph
from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.combinatorics.schubert_monomial_graph import SchubertMonomialGraph
from schubmult.symbolic import Expr, S, prod
from schubmult.utils._grid_print import GridPrint


class AntiRCGraph(SchubertMonomialGraph, GridPrint, CrystalGraph):
    _display_name = "AntiRCGraph"

    def __init__(self, rows_or_grid: Iterable[Iterable[int]] | np.ndarray, *, _is_copy: bool = False) -> None:
        if _is_copy:
            return
        self._perm = None
        self._word = None
        self._row_word_cache = None
        self._initialize_grid(rows_or_grid)

    def _initialize_grid(self, rows_or_grid: Iterable[Iterable[int]] | np.ndarray) -> None:
        if isinstance(rows_or_grid, np.ndarray):
            if rows_or_grid.ndim != 2:
                raise ValueError("AntiRCGraph grid must be a 2D ndarray")
            self._grid = np.array(rows_or_grid, dtype=np.uint8, copy=True)
            if np.any((self._grid != 0) & (self._grid != 1)):
                raise ValueError("AntiRCGraph grid entries must be 0 or 1")
            return

        rows = tuple(tuple(int(entry) for entry in row) for row in rows_or_grid)
        num_rows = len(rows)
        width = 0
        for row_index, row in enumerate(rows):
            row_label = num_rows - row_index
            if any(entry < row_label for entry in row):
                raise ValueError("Each anti reflection must be at least its anti row label")
            for entry in row:
                width = max(width, entry - row_label + 1)

        self._grid = np.zeros((num_rows, width), dtype=np.uint8)
        for row_index, row in enumerate(rows):
            row_label = num_rows - row_index
            for entry in row:
                self._grid[row_index, entry - row_label] = 1

    @property
    def rows(self) -> int:
        return self._grid.shape[0]

    @property
    def cols(self) -> int:
        return self._grid.shape[1]

    def __len__(self) -> int:
        return self.rows

    @cached_property
    def _row_words(self) -> tuple[tuple[int, ...], ...]:
        rows = []
        for row_index in range(self.rows):
            row_label = self.rows - row_index
            cols = np.flatnonzero(self._grid[row_index])
            rows.append(tuple(sorted((row_label + int(col) for col in cols), reverse=True)))
        return tuple(rows)

    def __iter__(self) -> Iterator[tuple[int, ...]]:
        return iter(self._row_words)

    def __getitem__(self, key: int | slice | tuple[int | slice, int | slice]):
        if isinstance(key, int):
            return self._row_words[key]

        if isinstance(key, tuple):
            i, j = key
            if isinstance(i, slice):
                return tuple(self[a, j] for a in range(self.rows)[i])
            if isinstance(j, slice):
                return tuple(self[i, b] for b in range(self.cols)[j])
            if not isinstance(i, int) or not isinstance(j, int):
                raise ValueError(f"Bad indexing {key=}")
            if i < 0 or j < 0 or i >= self.rows or j >= self.cols:
                return None
            if self._grid[i, j]:
                return self.rows - i + j
            return None

        if isinstance(key, slice):
            return self._row_words[key]

        raise ValueError(f"Bad indexing {key=}")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AntiRCGraph):
            return False
        return np.array_equal(self._grid, other._grid)

    def __hash__(self) -> int:
        return hash(self._grid.tobytes())

    def copy(self) -> AntiRCGraph:
        new_graph = AntiRCGraph(None, _is_copy=True)
        new_graph._grid = self._grid.copy()
        new_graph._perm = self._perm
        new_graph._word = self._word
        return new_graph

    def has_element(self, i: int, j: int) -> bool:
        if i <= 0 or j <= 0:
            return False
        if i > self.rows or j > self.cols:
            return False
        return bool(self._grid[i - 1, j - 1])

    @cached_property
    def anti_reduced_word(self) -> tuple[int, ...]:
        return tuple(entry for row in reversed(self._row_words) for entry in row)

    @cached_property
    def perm_word(self) -> tuple[int, ...]:
        # Read anti grid natively: right-to-left within each row, top-to-bottom across rows.
        word = []
        for i in range(self.rows):
            row_label = self.rows - i
            for j in range(self.cols - 1, -1, -1):
                if self._grid[i, j]:
                    word.append(row_label + j)
        return tuple(word)

    @property
    def reduced_word(self) -> tuple[int, ...]:
        return self.perm_word

    @property
    def anti_perm_word(self) -> tuple[int, ...]:
        return self.anti_reduced_word

    @cached_property
    def anti_compatible_sequence(self) -> tuple[int, ...]:
        sequence = []
        for row_label, row in enumerate(reversed(self._row_words), start=1):
            sequence.extend([row_label] * len(row))
        return tuple(sequence)

    @cached_property
    def anti_permutation(self) -> Permutation:
        return Permutation.ref_product(*self.anti_reduced_word)

    @property
    def perm(self) -> Permutation:
        return ~self.anti_permutation

    @property
    def reflection_view(self) -> RCGraph:
        return RCGraph(tuple(reversed(self._row_words)))

    def to_rc_graph(self) -> RCGraph:
        return self.reflection_view

    @classmethod
    def from_rc_graph(cls, rc: RCGraph) -> AntiRCGraph:
        return cls(tuple(reversed(tuple(rc))))

    @classmethod
    def from_reduced_anticompatible(
        cls,
        word: Sequence[int],
        seq: Sequence[int],
        length: int | None = None,
    ) -> AntiRCGraph:
        if len(word) != len(seq):
            raise ValueError("word and seq must have the same length")
        if any(entry <= 0 for entry in word) or any(row <= 0 for row in seq):
            raise ValueError("word and seq entries must be positive")
        if any(entry < row for entry, row in zip(word, seq, strict=False)):
            raise ValueError("anti reflections must be at least their anti row labels")

        inferred_length = max(seq, default=0)
        if length is None:
            length = inferred_length
        if length < inferred_length:
            raise ValueError("length must be at least max(seq)")

        rows = [[] for _ in range(length)]
        for entry, row_label in zip(word, seq, strict=False):
            rows[length - row_label].append(int(entry))
        canonical_rows = tuple(tuple(sorted(row, reverse=True)) for row in rows)
        return cls(canonical_rows)

    @property
    def anti_is_valid(self) -> bool:
        return all(entry >= row for entry, row in zip(self.anti_reduced_word, self.anti_compatible_sequence, strict=False)) and self.to_rc_graph().is_valid

    def as_reduced_anticompatible(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return self.anti_reduced_word, self.anti_compatible_sequence

    @cached_property
    def crystal_weight(self) -> tuple[int, ...]:
        return tuple(len(row) for row in self._row_words)

    @cached_property
    def length_vector(self) -> tuple[int, ...]:
        return self.crystal_weight

    def crystal_length(self) -> int:
        return self.rows

    def normalize(self) -> AntiRCGraph:
        return type(self).from_rc_graph(self.to_rc_graph().normalize())

    def polyvalue(self, x, y=None, **_kwargs) -> Expr:
        ret = S.One
        for i in range(self.rows):
            row_label = self.rows - i
            row_cols = np.flatnonzero(self._grid[i])
            if y is None:
                ret *= x[row_label] ** len(row_cols)
            else:
                # In anti coordinates, column j contributes y[j + 1] under row reflection.
                ret *= prod(x[row_label] - y[int(j) + 1] for j in row_cols)
        return ret

    def left_zero_act(self) -> set[AntiRCGraph]:
        return {type(self).from_rc_graph(rc) for rc in self.to_rc_graph().right_zero_act()}

    def right_zero_act(self) -> set[AntiRCGraph]:
        return self.left_zero_act()

    def antiaut(self) -> AntiRCGraph:
        return type(self)(self._grid[::-1, :].copy())

    def vertical_cut(self, row: int) -> tuple[AntiRCGraph, AntiRCGraph]:
        front_rc, back_rc = self.to_rc_graph().vertical_cut(self.rows - row)
        # Anti orientation swaps the cut factors relative to RC orientation.
        return type(self).from_rc_graph(back_rc), type(self).from_rc_graph(front_rc)

    def product(self, other: SchubertMonomialGraph) -> dict[AntiRCGraph, int]:
        other_rc = other.to_rc_graph() if isinstance(other, AntiRCGraph) else other
        product = other_rc.product(self.to_rc_graph())
        return {type(self).from_rc_graph(rc): coeff for rc, coeff in product.items()}

    def lowering_operator(self, row: int) -> AntiRCGraph | None:
        flipped = self.to_rc_graph().raising_operator(self.rows - row)
        if flipped is None:
            return None
        return type(self).from_rc_graph(flipped)

    def raising_operator(self, row: int) -> AntiRCGraph | None:
        flipped = self.to_rc_graph().lowering_operator(self.rows - row)
        if flipped is None:
            return None
        return type(self).from_rc_graph(flipped)

    @property
    def max_reflection(self) -> int:
        return max((entry for row in self for entry in row), default=0)

    def disjoint_union(self, anti_rc: AntiRCGraph) -> AntiRCGraph:
        if self.rows != anti_rc.rows:
            raise ValueError("Anti RC graphs must have the same number of rows")
        if self.perm.inv == 0:
            return anti_rc
        shift = self.max_reflection
        new_rows = self.rows + shift
        new_cols = max(self.cols, anti_rc.cols + shift)
        new_grid = np.zeros((new_rows, new_cols), dtype=np.uint8)

        # Growing the anti graph by `shift` rows corresponds to top-padding both inputs.
        new_grid[:self.rows, : self.cols] |= self._grid
        new_grid[:anti_rc.rows, shift : shift + anti_rc.cols] |= anti_rc._grid

        return type(self)(new_grid)


    def squash_product(self, anti_rc: AntiRCGraph) -> AntiRCGraph:
        dju = anti_rc.to_rc_graph().disjoint_union(self.to_rc_graph())
        return type(self).from_rc_graph(dju).vertical_cut(dju.rows - self.rows)[1]

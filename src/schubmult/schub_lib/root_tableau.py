import copy
from functools import cached_property
from typing import Any, Optional

import numpy as np
from sympy import pretty_print

from schubmult.utils._grid_print import GridPrint

from .crystal_graph import CrystalGraph
from .nilplactic import NilPlactic
from .perm_lib import Permutation
from .plactic import Plactic
from .rc_graph import RCGraph

# specific skew tableaux
# subword
# dual
# recti/fies to a specific subword
#subword subtableaux

# Dominant formula: number of standard tableaux of shape w/mu that recitfy to a subword tableau of shape u for a fixed tableau of shape w
# skew tableau behave dominant correctly
# w/mu subword tableau that are a specific standard tableaux for v

# we can do crazy crystal stuff

def _length_of_row(grid, row):
    return len([c for c in grid[row, :] if c is not None])

def _count_boxes(grid, spot=(0, 0)):
    count = 0
    for i in range(spot[0] ,grid.shape[0]):
        for j in range(spot[1] if i == spot[0] else 0, grid.shape[1]):
            cell = grid[i, j]
            if cell is not None:
                count += 1
    return count

def _word_from_grid(grid0, spot = (0, 0)):
    grid = copy.deepcopy(grid0)
    word = []
    def _recurse_grid():
        nonlocal grid, word
        max_r, max_c = -1, -1
        for i in range(grid.shape[0] - 1, spot[0] - 1, -1):
            L = _length_of_row(grid, i)
            if i > spot[0] or L >= spot[1]:
                cell = grid[i, L - 1] if L > 0 else None
                if cell is not None:
                    root = cell[0]
                    if root[1] == root[0] + 1:
                        max_r, max_c = i, L - 1
                        break
        if max_r == -1 and max_c == -1:
            pretty_print(grid)
            raise ValueError("No valid root found")
        root = grid[max_r, max_c][0]
        assert root[1] == root[0] + 1
        word.append(root[0])
        grid = _root_shift(root[0])(grid)
        grid[max_r, max_c] = None
        return _count_boxes(grid, spot=spot)
    while _recurse_grid() > 0:
        pass
    assert _count_boxes(grid) == 0, f"Grid should be empty after extraction, but found {_count_boxes(grid)} boxes, {grid=}."
    assert len(word) == _count_boxes(grid0)
    return tuple(reversed(word))


def _root_shift(root):
    """
    Return a callable shift(grid_slice) -> new_grid_slice that applies the
    appropriate root-reflection action to every non-None cell of the input
    object-array slice. Uses Permutation.ref_product(...) as the reflection
    provider and is defensive about act_root signatures.
    """
    # choose a permutation-reflector if we can
    sref = None
    if isinstance(root, int):
        sref = Permutation.ref_product(int(root))
    elif isinstance(root, (tuple, list)) and len(root) >= 1:
        # use first entry as the simple-reflection index (safe fallback)
        sref = Permutation.reflection(root)

    def _shift(grid_slice):
        # expect a numpy object-array (or something indexable)
        grid_slice = np.asarray(grid_slice, dtype=object)
        out = np.empty(grid_slice.shape, dtype=object)

        # support 1-D and 2-D slices (ndindex yields tuples whose length equals ndim)
        if grid_slice.ndim == 1:
            for i0 in range(grid_slice.shape[0]):
                cell = grid_slice[i0]
                if cell is None:
                    out[i0] = None
                    continue
                root_cell, letter = cell
                new_root = root_cell
                if sref is not None:
                    new_root = sref.act_root(root_cell)
                new_root = tuple(int(x) for x in new_root)
                out[i0] = (new_root, letter)
        else:
            for i0, j0 in np.ndindex(grid_slice.shape):
                cell = grid_slice[i0, j0]
                if cell is None:
                    out[i0, j0] = None
                    continue
                root_cell, letter = cell
                new_root = root_cell
                if sref is not None:
                    new_root = sref.act_root(*root_cell)
                new_root = tuple(int(x) for x in new_root)
                out[i0, j0] = (new_root, letter)
        return out

    return _shift


class RootTableau(CrystalGraph, GridPrint):
    """
    Root tableau with dual knuth equivalence
    """

    def __hash__(self) -> int:
        return hash(self._hasher)


    @classmethod
    def root_insert_rsk(cls, reduced_word, compatible_seq):
        _perm = Permutation.ref_product(*reduced_word)
        word, word2 = (), ()
        spunkle = len(_perm)
        w0 = Permutation.w0(spunkle)
        rev_word = [w0[a-1] for a in reduced_word]

        for idx, letter in enumerate(rev_word):
            letter2 = idx + 1
            word, word2 = NilPlactic._ed_insert_rsk(word, word2, int(letter), int(letter2) if letter2 is not None else None)
        num_rows = len(word2)
        num_cols = max(len(r) for r in word2)
        grid = np.empty((num_rows, num_cols), dtype=object)
        for r in range(num_rows):
            for c in range(num_cols):
                if c < len(word2[r]):
                    # store pair (root_tuple, letter) as before
                    grid[r, c] = (_perm.right_root_at(word2[r][c] - 1, word=reduced_word), compatible_seq[word2[r][c] - 1])
                else:
                    grid[r, c] = None

        return cls(grid)

    # skew tableaux are subword
    @classmethod
    def from_rc_graph(cls, rc: RCGraph):
        reduced_word = rc.perm_word
        compatible_seq = []
        for i in range(len(rc)):
            compatible_seq.extend([i+1] * len(rc[i]))
        return cls.root_insert_rsk(reduced_word, compatible_seq)

    def _index_of_box(self, row, col):
        return sum(self.length_of_row(r) for r in range(row)) + col


    def delete_box(self, i, j):
        if self[i, j] is None:
            raise ValueError(f"Missing box at position {(i, j)} cannot be deleted")

        # use a deep copy so we never mutate the original tableau's objects/views
        
        def _delete_box_jdt():
            """
            Delete the box at (i,j) (0-indexed) by creating a hole, applying the
            appropriate root-shift to the regions to the left/above the hole, then
            performing the downward/rightward jeu-de-taquin slide (push into hole).

            Returns a new numpy object-array (deep-copied) with the box removed and
            shape trimmed of trailing empty rows/columns.
            """
            nonlocal i, j
            new_grid = copy.deepcopy(self._root_grid)

            rows, cols = new_grid.shape
            if new_grid[i, j] is None:
                raise ValueError(f"No box at position {(i,j)} to delete")

            # root = new_grid[i, j][0]
            new_grid[i, j] = None

            # apply root-shift to region above the deleted box (all columns) and to
            # the part of the same row left of the hole; be defensive about shapes
            #new_grid = _root_shift_(root)(new_grid)
            new_grid[:i + 1, :j +1] = _root_shift(self.letter_at(i, j))(new_grid[:i + 1, :j + 1])
            # perform down/right jeu-de-taquin (push boxes into the hole)
            rows, cols = new_grid.shape
            while True:
                down = new_grid[i + 1, j] if (i + 1 < rows and j < cols) else None
                right = new_grid[i, j + 1] if (j + 1 < cols) else None

                if down is None and right is None:
                    break

                if down is None and right is not None:
                    new_grid[i, j] = right
                    j += 1
                elif right is None and down is not None:
                    new_grid[i, j] = down
                    i += 1
                else:
                    # pick smaller by recording letter (second component)
                    down_val = down[1] if (isinstance(down, tuple) and len(down) > 1) else down
                    right_val = right[1] if (isinstance(right, tuple) and len(right) > 1) else right
                    if down_val <= right_val:
                        new_grid[i, j] = down
                        i += 1
                    else:
                        new_grid[i, j] = right
                        j += 1

            # remove the final moved box (the hole reaches an outer cell)
            new_grid[i, j] = None

            # trim trailing empty columns
            def _col_is_empty(g, c):
                return all(g[r, c] is None for r in range(g.shape[0]))

            def _row_is_empty(g, r):
                return all(g[r, c] is None for c in range(g.shape[1]))

            # trim rightmost empty columns
            while new_grid.shape[1] > 0 and _col_is_empty(new_grid, new_grid.shape[1] - 1):
                # build smaller array without last column
                rcount, ccount = new_grid.shape
                if ccount == 1:
                    new_grid = np.empty((rcount, 0), dtype=object)
                    break
                tmp = np.empty((rcount, ccount - 1), dtype=object)
                for rr in range(rcount):
                    for cc in range(ccount - 1):
                        tmp[rr, cc] = new_grid[rr, cc]
                new_grid = tmp

            # trim bottom empty rows
            while new_grid.shape[0] > 0 and _row_is_empty(new_grid, new_grid.shape[0] - 1):
                rcount, ccount = new_grid.shape
                if rcount == 1:
                    new_grid = np.empty((0, ccount), dtype=object)
                    break
                tmp = np.empty((rcount - 1, ccount), dtype=object)
                for rr in range(rcount - 1):
                    for cc in range(ccount):
                        tmp[rr, cc] = new_grid[rr, cc]
                new_grid = tmp

            return new_grid

        return self.__class__(_delete_box_jdt())


    def __getitem__(self, key: Any) -> Any:
        return self._root_grid[key]

    @cached_property
    def rows(self):
        return self._root_grid.shape[0]

    @cached_property
    def cols(self):
        return self._root_grid.shape[1]

    # @property
    # def compatible_sequence(self):
    #     return self._plactic.reverse_rsk(self._index_tableau)

    @property
    def is_valid(self):
        return self.rc_graph.is_valid

    @property
    def reduced_word(self):
        return _word_from_grid(self._root_grid)

    def letter_at(self, row, col):
        return _word_from_grid(self._root_grid, spot=(row, col))[0]
    # @property
    # def reduced_word(self):
    #     return self._red_plactic.reverse_rsk(self._index_tableau)

    def __init__(self, grid=None):
        self._root_grid = copy.deepcopy(grid) if grid is not None else np.empty((0, 0), dtype=object)
        self._hasher = tuple(tuple(tuple(b) for b in a if b is not None) for a in self._root_grid if a is not None)


    @property
    def weight_tableau(self):
        _word = tuple([a for a in row if a is not None] for row in self._root_grid)
        return Plactic(_word)

    def epsilon(self, index):
        return self._weight_tableau.epsilon(index)

    def raising_operator(self, index):
        new_plactic = self._weight_tableau.raising_operator(index)
        if new_plactic is None:
            return None
        return RootTableau(self._base_word, new_plactic)

    def lowering_operator(self, index):
        new_plactic = self._weight_tableau.lowering_operator(index)
        if new_plactic is None:
            return None
        new_tab = RootTableau(self._base_word, new_plactic)
        if new_tab.rc_graph is None or not new_tab.rc_graph.is_valid:
            return None
        return new_tab

    @property
    def rc_graph(self):
        rows = []
        for i in range(len(self._base_word)):
            if i == 0 or (i > 0 and self._base_word[i] > self._base_word[i - 1]):
                rows.append([])
            rows[-1].append(self._base_word[i])
        raise_seq = self._weight_tableau.to_highest_weight(length=self.crystal_length())[1]
        try:
            rc = RCGraph([tuple(row) for row in rows]).normalize()
            return rc.reverse_raise_seq(raise_seq)
        except ValueError:
            return None

    def crystal_length(self):
        return len(self._perm.trimcode)

    @property
    def crystal_weight(self):
        return self._plactic.crystal_weight()

    # def reduced_word(self):
    #     recording_tableau = self.index_tableau
    #     permutation_of_roots = self.reverse_rsk(recording_tableau)
    #     roots = list(reversed([self._perm.right_root_at(i - 1) for i in permutation_of_roots]))
    #     word = []
    #     for i in range(len(roots)):
    #         word.append(roots[0][0])
    #         sref = Permutation.ref_product(roots[0][0])
    #         roots = [sref.act_root(*r) for r in roots[1:]]
    #     return tuple(reversed(word))




    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RootTableau):
            return False
        return self._hasher == other._hasher
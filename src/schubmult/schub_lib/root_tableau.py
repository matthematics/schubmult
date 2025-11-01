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
# subword subtableaux

# Dominant formula: number of standard tableaux of shape w/mu that recitfy to a subword tableau of shape u for a fixed tableau of shape w
# skew tableau behave dominant correctly
# w/mu subword tableau that are a specific standard tableaux for v

# we can do crazy crystal stuff


def _length_of_row(grid, row):
    return len([c for c in grid[row, :] if c is not None])


def _count_boxes(grid):
    count = 0
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            cell = grid[i, j]
            if cell is not None:
                count += 1
    return count


def _root_compare(root1, root2):
    if root1 == root2:
        return 2
    if root1[1] == root2[1] and root1[0] != root2[0]:
        return 1
    if root1[0] == root2[0] and root1[1] != root2[1]:
        return 1
    if root1[0] == root2[1] or root1[1] == root2[0]:
        return -1
    return 0


def _word_from_grid(grid0, as_grid: Optional[bool] = False, as_ordering: Optional[bool] = False, with_compatible_seq: Optional[bool] = False) -> Any:
    """
    Two modes:
      - as_grid=True: return an object-array the same shape as grid0 where each
        occupied cell contains the recording letter (cell[1]) and empty cells
        are None.
      - as_grid=False: reconstruct the reduced word (sequence of simple-reflection
        indices) by repeatedly:
          * finding the occupied cell whose recording letter (cell[1]) is maximal
            and, among those, is farthest to the right (largest column index;
            break ties by largest row index),
          * popping that box, appending cell[0][0] to the collected word,
          * applying the root-shift corresponding to cell[0] to the region
            above the popped box and to the part of the same row left of the box,
          * repeating until no boxes remain.
        Returns a tuple(reversed(collected_letters)) to match the original reduced
        word orientation used elsewhere.
    """
    ARBITRARY_BIG_NUMBER = 1000
    def _flip_grid(grid):
        nonlocal ARBITRARY_BIG_NUMBER
        for ii in range(grid.shape[0]):
            for jj in range(grid.shape[1]):
                cell = grid[ii, jj]
                if cell is None:
                    continue
                root_cell, letter = cell
                new_root = (ARBITRARY_BIG_NUMBER - root_cell[1], ARBITRARY_BIG_NUMBER - root_cell[0])
                grid[ii, jj] = (new_root, letter)
    def _flip(i):
        return ARBITRARY_BIG_NUMBER - i
    if as_grid:
        index_val = _count_boxes(grid0)
        # Build an output array the same shape as the grid and place, at the
        # location where boxes are popped during the reconstruction procedure,
        # the letter that is "found" at that pop (we use the pop-rule below:
        # pick max recording value farthest right; the popped letter is the
        # first component of the root cell).
        grid = copy.deepcopy(np.asarray(grid0, dtype=object))
        _flip_grid(grid)
        out = np.full(grid.shape, None, dtype=object)
        ordering = np.full(grid.shape, None, dtype=object)

        def boxes_remaining_local(g):
            for ii in range(g.shape[0]):
                for jj in range(g.shape[1]):
                    if g[ii, jj] is not None:
                        return True
            return False

        while boxes_remaining_local(grid):
            # find maximal recording letter value
            max_val = None
            for ii in range(grid.shape[0]):
                for jj in range(grid.shape[1]):
                    cell = grid[ii, jj]
                    if cell is None:
                        continue
                    val = cell[1]
                    if max_val is None or val > max_val:
                        max_val = val
            if max_val is None:
                break

            # choose the rightmost (then bottom-most) cell with that recording value
            chosen = None
            for ii in range(grid.shape[0]):
                for jj in range(grid.shape[1]):
                    cell = grid[ii, jj]
                    if cell is None or cell[1] != max_val:
                        continue
                    if chosen is None:
                        chosen = (ii, jj)
                    else:
                        ci, cj = chosen
                        if jj > cj or (jj == cj and ii > ci):
                            chosen = (ii, jj)
            if chosen is None:
                break

            i, j = chosen
            cell = grid[i, j]
            root_cell, _letter = cell
            # place the popped letter into the output at the popped location
            out[i, j] = _flip(root_cell[1])

            ordering[i, j] = index_val - 1
            index_val -= 1
            # remove box and apply shift to left / above regions
            rd = int(root_cell[0]) if isinstance(root_cell, (tuple, list)) else int(root_cell)
            grid[i, j] = None
            grid = _root_shift(rd)(grid)

        return out if not as_ordering else ordering

    # reconstruct reduced word by repeated deletion
    grid = copy.deepcopy(grid0)
    _flip_grid(grid)
    word = []
    compatible_seq = []

    def boxes_remaining(g):
        for ii in range(g.shape[0]):
            for jj in range(g.shape[1]):
                if g[ii, jj] is not None:
                    return True
        return False

    while boxes_remaining(grid):
        # find maximal recording letter value
        max_val = None
        for ii in range(grid.shape[0]):
            for jj in range(grid.shape[1]):
                cell = grid[ii, jj]
                if cell is None:
                    continue
                val = cell[1]
                if max_val is None or val > max_val:
                    max_val = val

        if max_val is None:
            break  # no boxes

        # among cells with recording == max_val choose the one farthest to the right
        # (largest column index). Break ties by largest row index.
        chosen = None  # (i,j)
        for ii in range(grid.shape[0]):
            if chosen is not None:
                break
            for jj in range(grid.shape[1]):
                cell = grid[ii, jj]
                if cell is None:
                    continue
                if cell[1] != max_val:
                    continue
                if chosen is None:
                    chosen = (ii, jj)
                else:
                    ci, cj = chosen
                    # prefer larger column, then larger row
                    if jj == grid.shape[1] - 1 or grid[ii, jj + 1] is None:
                        chosen = (ii, jj)
                        break

        if chosen is None:
            # nothing found (shouldn't happen)
            break

        i, j = chosen
        cell = grid[i, j]
        
        root_cell, _letter = cell
        # append the first component of the root as the letter for the reduced word
        word.append(_flip(root_cell[1]))
        compatible_seq.append(_letter)

        # remove the box and reflect the remaining grid before continuing
        # apply root-shift corresponding to the popped root to (row[:j]) and ([:i, :])
        # Use the integer index of the root (first entry) as the shift key
        rd = int(root_cell[0])
        # set removed box to None first (so _root_shift sees correct shape)
        grid[i, j] = None
        grid = _root_shift(rd)(grid)

    # the algorithm collected letters in pop order; return reversed to match original orientation
    if with_compatible_seq:
        return tuple(reversed(word)), tuple(reversed(compatible_seq))
    return tuple(reversed(word))


def _root_shift(root, spots=None):
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
                if spots is not None and i0 not in spots:
                    out[i0] = grid_slice[i0]
                    continue
                cell = grid_slice[i0]
                if cell is None:
                    out[i0] = None
                    continue
                root_cell, letter = cell
                new_root = root_cell
                if sref is not None:
                    new_root = sref.act_root(*root_cell)
                    if new_root[0] > new_root[1]:
                        new_root = root_cell
                new_root = tuple(int(x) for x in new_root)
                out[i0] = (new_root, letter)
        else:
            for i0, j0 in np.ndindex(grid_slice.shape):
                if spots is not None and (i0, j0) not in spots:
                    out[i0, j0] = grid_slice[i0, j0]
                    continue
                cell = grid_slice[i0, j0]
                if cell is None:
                    out[i0, j0] = None
                    continue
                root_cell, letter = cell
                new_root = root_cell
                if sref is not None:
                    new_root = sref.act_root(*root_cell)
                    if new_root[0] > new_root[1]:
                        new_root = root_cell
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
        rev_word = [w0[a - 1] for a in reduced_word]

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
            compatible_seq.extend([i + 1] * len(rc[i]))
        return cls.root_insert_rsk(reduced_word, compatible_seq)

    def _index_of_box(self, row, col):
        return sum(self.length_of_row(r) for r in range(row)) + col

    def roots_before(self, row, col):
        order_grid = _word_from_grid(self._root_grid, as_ordering=True, as_grid=True)
        return [(i, j) for (i, j) in np.ndindex(order_grid.shape) if order_grid[i, j] is not None and order_grid[i, j] < order_grid[row, col]]

    def up_jdt_slide(self, row, col):
        if self[row, col] != None:
            raise ValueError("Can only slide from empty box")
        if self[row - 1, col] is None and self[row, col - 1] is None:
            raise ValueError("No boxes to slide from")
        new_grid = copy.deepcopy(self._root_grid)

        def _recurse():
            nonlocal row, col, new_grid
            if row == 0 or (col > 0 and row > 0 and new_grid[row - 1, col] is None):
                # slide from left
                new_grid[row, col] = new_grid[row, col - 1]
                new_grid[row, col - 1] = None
                col -= 1
            elif col == 0 or (col > 0 and row > 0 and new_grid[row, col - 1] is None):
                # slide from above
                new_grid[row, col] = new_grid[row - 1, col]
                new_grid[row - 1, col] = None
                row -= 1
            else:
                # both available, pick larger root
                root_above = new_grid[row - 1, col][1]
                root_left = new_grid[row, col - 1][1]
                if root_above >= root_left:
                    # above is larger or incomparable
                    new_grid[row, col] = new_grid[row - 1, col]
                    new_grid[row - 1, col] = None
                    row -= 1
                else:
                    # left is larger
                    new_grid[row, col] = new_grid[row, col - 1]
                    new_grid[row, col - 1] = None
                    col -= 1
            if row > 0 or col > 0:
                _recurse()

        _recurse()
        return RootTableau(new_grid)

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
        return _word_from_grid(self._root_grid, as_grid=False)

    @property
    def word_grid(self):
        return _word_from_grid(self._root_grid, as_grid=True)

    def letter_at(self, row, col):
        return self.word_grid[row, col]

    # @property
    # def reduced_word(self):
    #     return self._red_plactic.reverse_rsk(self._index_tableau)

    @cached_property
    def print_element(self):
        _printing_grid = copy.deepcopy(self._root_grid)
        for i in range(_printing_grid.shape[0]):
            for j in range(_printing_grid.shape[1]):
                cell = _printing_grid[i, j]
                if cell is None:
                    _printing_grid[i, j] = " "
        return RootTableau(_printing_grid)

    def __init__(self, grid):
        self._root_grid = copy.deepcopy(grid)
        self._hasher = tuple(tuple(tuple(b) for b in a if b is not None) for a in self._root_grid if a is not None)

    @property
    def weight_tableau(self):
        _word = tuple([a for a in row if a is not None] for row in self._root_grid)
        return Plactic(_word)

    def epsilon(self, index):
        return self._weight_tableau.epsilon(index)

    def raising_operator(self, index):
        # up = self.rc_graph.raising_operator(index)
        # if up is None:
        #     return None
        # new_grid = copy.deepcopy(self._root_grid)
        # root_map = dict([(up.left_to_right_inversions(i), up.left_to_right_inversion_coords(i)[0]) for i in range(up.perm.inv)])
        # for i, j in new_grid.
        raise NotImplementedError("RCGraph API not implemented will need to do it for real")
        up = self.rc_graph.raising_operator(index)
        if up is None:
            return None

        # deep copy so we don't mutate self
        new_grid = copy.deepcopy(self._root_grid)

        # build mapping from whatever keys the RCGraph produces to the new recording value
        # root_map keys/values depend on RCGraph API; be defensive when applying.
        root_map = {up.left_to_right_inversion(i): up.left_to_right_inversion_coords(i)[0] for i in range(up.perm.inv)}

        # iterate over every cell and replace the second component (recording letter)
        # according to root_map when possible. Tuples are immutable so we construct a new tuple.
        
        for ii, jj in np.ndindex(new_grid.shape):
            cell = new_grid[ii, jj]
            if cell is None:
                continue
            root_cell, _ = cell
            try:
                new_letter = root_map[root_cell]
            except Exception:
                raise ValueError(f"Error looking up {root_cell} in {root_map} from RC graph {up} {up.perm=} {self=}")

            new_grid[ii, jj] = (root_cell, new_letter)

        return RootTableau(new_grid)

    def lowering_operator(self, index):

        raise NotImplementedError("RCGraph API not implemented will need to do it for real")
        down = self.rc_graph.lowering_operator(index)
        if down is None:
            return None
        # deep copy so we don't mutate self
        new_grid = copy.deepcopy(self._root_grid)

        # build mapping from whatever keys the RCGraph produces to the new recording value
        # root_map keys/values depend on RCGraph API; be defensive when applying.
        root_map = {down.left_to_right_inversion(i): down.left_to_right_inversion_coords(i)[0] for i in range(down.perm.inv)}

        # iterate over every cell and replace the second component (recording letter)
        # according to root_map when possible. Tuples are immutable so we construct a new tuple.
        for ii, jj in np.ndindex(new_grid.shape):
            cell = new_grid[ii, jj]
            if cell is None:
                continue
            root_cell, _ = cell
            # try several sensible lookups (be permissive about key types)
            new_letter = root_map[root_cell]

            new_grid[ii, jj] = (root_cell, new_letter)

        return RootTableau(new_grid)

    @property
    def rc_graph(self):
        reduced_word, compatible_seq = _word_from_grid(self._root_grid, with_compatible_seq=True)
        rows = []
        assert len(reduced_word) == len(compatible_seq)
        print([(a,b) for a,b in zip(reduced_word, compatible_seq)])
        for i, a in enumerate(compatible_seq):
            if a > len(rows):
                while len(rows) < a:
                    rows.append(())
            rows[a-1] = (*rows[a-1], reduced_word[i])
        return RCGraph(tuple(rows)).normalize()

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

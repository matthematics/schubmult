import copy
import logging
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
def _plactic_raising_operator(word, i):
    word = [*word]
    opening_stack = []
    closing_stack = []
    for index in range(len(word)):
        if word[index] == i + 1:
            opening_stack.append(index)
        elif word[index] == i:
            if len(opening_stack) > 0:
                opening_stack.pop()
            else:
                closing_stack.append(index)
    if len(opening_stack) == 0:
        return None
    index_to_change = opening_stack[0]
    word[index_to_change] = i
    return tuple(word)

def _is_valid_outer_corner(grid: np.ndarray, i: int, j: int) -> bool:
    """
    Outer-corner predicate used by up_jdt_slide.
    Accepts hole positions that may extend the grid (i==rows or j==cols).
    A position is valid if:
      - there is a box above AND a box to the left, OR
      - at least one of those exists and the hole is on/extends the outer boundary.
    """
    rows, cols = grid.shape
    # treat positions outside current array as empty slots (they must be extended before sliding)
    if 0 <= i < rows and 0 <= j < cols and grid[i, j] is not None:
        return False
    if i >= rows or j >= cols:
        return False
    up = grid[i - 1, j] if i - 1 >= 0 else "bob" if i == 0 else None
    left = grid[i, j - 1] if j - 1 >= 0 else "bing" if j == 0 else None
    return grid[i, j] is None and (up is not None and left is not None) and not (i == 0 and j == 0)
    # consider hole on or beyond boundary as "outer


def _is_valid_inner_corner(grid: np.ndarray, i: int, j: int) -> bool:
    """
    Inner-corner predicate used by down_jdt_slide.
    Valid when the hole is inside the grid (must not extend grid) and there
    is a box below or to the right with the same boundary rules mirrored.
    """
    rows, cols = grid.shape
    if grid[i, j] is not None:
        return False
    if i < 0 or j < 0:
        return False
    if i >= rows or j >= cols:
        return False
    down = grid[i + 1, j] if i + 1 < rows else "bob" if i == rows - 1 else None
    right = grid[i, j + 1] if j + 1 < cols else "bing" if j == cols - 1 else None
    return grid[i, j] is None and (down is not None and right is not None) and not (i == rows - 1 and j == cols - 1)


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


_logger = logging.getLogger(__name__)
if not _logger.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    _logger.addHandler(_h)
_logger.setLevel(logging.DEBUG)


def _validate_grid(grid: np.ndarray) -> None:
    """Lightweight validation of a root-grid; raises on malformed cells."""
    if grid is None:
        raise RuntimeError("_validate_grid: grid is None")
    if not hasattr(grid, "shape"):
        raise RuntimeError("_validate_grid: grid missing shape")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            cell = grid[i, j]
            if cell is None:
                continue
            if not isinstance(cell, (tuple, list)) or len(cell) < 2:
                raise RuntimeError(f"_validate_grid: bad cell at {(i, j)}: {repr(cell)}")
            root, letter = cell[0], cell[1]
            try:
                int(letter)
            except Exception:
                raise RuntimeError(f"_validate_grid: non-int letter at {(i, j)}: {repr(letter)}")
            if not (isinstance(root, (tuple, list, int))):
                raise RuntimeError(f"_validate_grid: unexpected root at {(i, j)}: {repr(root)}")


def _snap_grid(grid: np.ndarray):
    """Return a compact, JSON-like snapshot of the grid for debug messages."""
    try:
        return [[None if grid[i, j] is None else (grid[i, j][0], int(grid[i, j][1])) for j in range(grid.shape[1])] for i in range(grid.shape[0])]
    except Exception:
        return repr(grid)


def _root_map(rc1, rc2):
    # takes roots of rc1 return roots of rc2 (dct)

    rw1 = rc1
    rw2 = rc2
    perm = Permutation.ref_product(*rw1)
    return {perm.right_root_at(i, word=rw1): perm.right_root_at(i, word=rw2) for i in range(perm.inv)}


class RootTableau(CrystalGraph, GridPrint):
    """
    Root tableau with dual knuth equivalence
    """

    def __hash__(self) -> int:
        return hash(self._hasher)

    @property
    def perm(self):
        return Permutation.ref_product(*self.reduced_word)

    # preserved by the crystal operators
    @property
    def edelman_greene_invariant(self):
        w0 = Permutation.w0(len(self.perm))
        rev_word = [w0[r - 1] for r in self.reduced_word]
        return NilPlactic().ed_insert(*rev_word)

    @property
    def shape(self):
        return tuple(_length_of_row(self._root_grid, r) for r in range(self._root_grid.shape[0]) if _length_of_row(self._root_grid, r) > 0)

    @classmethod
    def root_insert_rsk(cls, reduced_word, compatible_seq):
        _perm = Permutation.ref_product(*reduced_word)
        #word, word2 = (), ()
        spunkle = len(_perm)
        w0 = Permutation.w0(spunkle)
        rev_word = [w0[a - 1] for a in reduced_word]
        NP, P = NilPlactic.ed_insert_rsk(*rev_word)
        num_rows = len(NP._word)
        try:
            num_cols = len(NP._word[0])
        except IndexError:
            num_cols = 0
        grid = np.empty((num_rows, num_cols), dtype=object)
        for box in P.iter_boxes:
            grid[box] = (_perm.right_root_at(P[box] - 1, word=reduced_word), compatible_seq[P[box] - 1])
 
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

    @property
    def perm(self):
        return Permutation.ref_product(*self.reduced_word)

    def rectify(self, randomized=False):
        import random

        cur = self
        while True:
            inner_corners = tuple(cur.iter_inner_corners)
            if not inner_corners:
                break
            if randomized:
                cur = cur.down_jdt_slide(*random.choice(inner_corners))
            else:
                cur = cur.down_jdt_slide(*next(iter(inner_corners)))
        return cur

    def up_jdt_slide(self, row, col, check=False):
        if not _is_valid_outer_corner(self._root_grid, row, col):
            raise ValueError("Can only slide from valid outer corner")
        new_grid = copy.deepcopy(self._root_grid)

        def _recurse():
            nonlocal row, col, new_grid
            # _logger.debug(f"{row=} {col=}")
            #  pretty_print(RootTableau(new_grid))
            left = new_grid[row, col - 1] if col > 0 else None
            # slide from left
            up = new_grid[row - 1, col] if row > 0 else None
            # slide from above
            if up is None and left is None:
                return
            if left is None:
                new_grid[row, col] = up
                new_grid[row - 1, col] = None
                row -= 1
            elif up is None:
                new_grid[row, col] = left
                new_grid[row, col - 1] = None
                col -= 1
            else:
                # both available, pick larger root
                root_above = new_grid[row - 1, col][1]
                root_left = new_grid[row, col - 1][1]
                if root_above >= root_left:
                    # above is larger or incomparable
                    new_grid[row, col] = up
                    new_grid[row - 1, col] = None
                    row -= 1
                else:
                    # left is larger
                    new_grid[row, col] = left
                    new_grid[row, col - 1] = None
                    col -= 1
            _recurse()

        _recurse()
        new_grid[row, col] = None
        ret = RootTableau(new_grid)
        if check:
            assert ret.rc_graph == self.rc_graph, "up_jdt_slide does not preserve RC graph"
            assert ret.weight_tableau == self.weight_tableau, "up_jdt_slide does not preserve tableau shape"
        assert self.edelman_greene_invariant == ret.edelman_greene_invariant
        return ret

    def down_jdt_slide(self, row, col, check=False):
        """
        Perform a downward/rightward jeu-de-taquin slide starting from the given
        (row, col) hole (0-indexed). Boxes from below or to the right are moved
        into the hole, preferring the smaller recording letter when both exist.
        Returns a new RootTableau (does not mutate self).
        """
        if not _is_valid_inner_corner(self._root_grid, row, col):
            raise ValueError("Can only slide from valid inner corner")

        new_grid = copy.deepcopy(self._root_grid)

        def _recurse():
            nonlocal row, col, new_grid
            # _logger.debug(f"{row=} {col=}")
            #  pretty_print(RootTableau(new_grid))
            right = new_grid[row, col + 1] if col < self.cols - 1 else None
            # slide from left
            down = new_grid[row + 1, col] if row < self.rows - 1 else None
            # slide from above
            if right is None and down is None:
                return
            if right is None:
                new_grid[row, col] = down
                new_grid[row + 1, col] = None
                row += 1
            elif down is None:
                new_grid[row, col] = right
                new_grid[row, col + 1] = None
                col += 1
            else:
                # both available, pick larger root
                root_below = new_grid[row + 1, col][1]
                root_right = new_grid[row, col + 1][1]
                if root_below <= root_right:
                    # below is larger or incomparable
                    new_grid[row, col] = down
                    new_grid[row + 1, col] = None
                    row += 1
                else:
                    # right is larger
                    new_grid[row, col] = right
                    new_grid[row, col + 1] = None
                    col += 1
            _recurse()

        _recurse()
        new_grid[row, col] = None

        ret = RootTableau(new_grid)
        if check:
            assert ret.rc_graph == self.rc_graph, "down_jdt_slide does not preserve RC graph"
            assert ret.weight_tableau == self.weight_tableau, "down_jdt_slide does not preserve weight tableau"
        assert self.edelman_greene_invariant == ret.edelman_greene_invariant
        return ret

    def __getitem__(self, key: Any) -> Any:
        return self._root_grid[key]

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, RootTableau):
            return False
        return self._hasher == other._hasher

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
    def iter_boxes(self):
        for i in range(self.rows):
            for j in range(self.cols):
                cell = self._root_grid[i, j]
                if cell is not None:
                    yield (i, j)

    @property
    def iter_outer_corners(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if _is_valid_outer_corner(self._root_grid, i, j):
                    yield (i, j)

    @property
    def iter_inner_corners(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if _is_valid_inner_corner(self._root_grid, i, j):
                    yield (i, j)

    @property
    def is_valid(self):
        return self.rc_graph.is_valid

    @property
    def reduced_word(self):
        return _word_from_grid(self._root_grid, as_grid=False)

    @property
    def compatible_sequence(self):
        return _word_from_grid(self._root_grid, with_compatible_seq=True)[1]

    @property
    def word_grid(self):
        return _word_from_grid(self._root_grid, as_grid=True)

    @property
    def grid_word(self):
        return tuple([self.word_grid[box] for box in self.iter_boxes_row_word_order()])

    @property
    def order_grid(self):
        return _word_from_grid(self._root_grid, as_grid=True, as_ordering=True)

    def letter_at(self, row, col):
        return self.word_grid[row, col]

    # @property
    # def reduced_word(self):
    #     return self._red_plactic.reverse_rsk(self._index_tableau)

    # @cached_property
    # def print_element(self):
    #     _printing_grid = copy.deepcopy(self._root_grid)
    #     for i in range(_printing_grid.shape[0]):
    #         for j in range(_printing_grid.shape[1]):
    #             cell = _printing_grid[i, j]
    #             if cell is None:
    #                 _printing_grid[i, j] = " "
    #     return RootTableau(_printing_grid, print_only=True)

    def __init__(self, grid, print_only=False):
        self._root_grid = copy.deepcopy(grid)
        self._hasher = tuple(tuple(tuple(b) for b in a if b is not None) for a in self._root_grid if a is not None)
        if not print_only:
            for ind in self.iter_boxes:
                if self[ind][0] is not None:
                    assert self[ind] == (self.perm.right_root_at(self.order_grid[ind], word=self.reduced_word), self.compatible_sequence[self.order_grid[ind]]), (
                        f"RootTableau init: inconsistent root at {ind}: {self[ind][0]=} vs {self.perm.right_root_at(self.order_grid[ind], word=self.reduced_word)=} {self.reduced_word=}"
                    )
                    # assert self.perm == Permutation.ref_product(*self.grid_word), f"{self.reduced_word=} {self.grid_word=}"
                    # for index, box in enumerate(reversed(list(self.iter_boxes_row_word_order()))):
                    #     assert self[box][0] == self.perm.right_root_at(index, word=list(reversed(self.grid_word)))

    @property
    def weight_tableau(self):
        # rows = []
        # for row in self._root_grid:
        #     rows.append([])
        #     for cell in row:
        #         if cell is None:
        #             rows[-1].append(0)
        #         else:
        #             rows[-1].append(cell[1])

        return Plactic().rs_insert(*self.row_word)

    def epsilon(self, index):
        return self._weight_tableau.epsilon(index)

    @property
    def row_word(self):
        word = []
        for r in range(self.rows - 1, -1, -1):
            for c in range(self.cols):
                cell = self._root_grid[r, c]
                if cell is not None:
                    word.append(cell[1])
        return tuple(word)

    @property
    def root_row_word(self):
        word = []
        for r in range(self.rows - 1, -1, -1):
            for c in range(self.cols):
                cell = self._root_grid[r, c]
                if cell is not None:
                    root_cell, _letter = cell
                    word.append(root_cell)
        return tuple(word)

    def right_root_at(self, i):
        for ind in self.iter_boxes:
            if self.order_grid[ind] is not None and self.order_grid[ind] == i:
                return self.perm.right_root_at(self.order_grid[ind], word=self.reduced_word)
        return None

    def iter_boxes_row_word_order(self):
        for i in range(self._root_grid.shape[0] - 1, -1, -1):
            for j in range(self._root_grid.shape[1]):
                if self[i, j] is not None:
                    yield (i, j)

    def raising_operator(self, i):
        """Crystal raising operator e_i on the root tableau"""
        row = i
        rc = self.rc_graph
        if row >= len(rc):
            return None
        row_i = [*rc[row - 1]]
        row_ip1 = [*rc[row]]

        # pair the letters
        pairings = []
        unpaired = []
        unpaired_b = [*row_ip1]

        for letter in row_i:
            st = [letter2 for letter2 in unpaired_b if letter2 > letter]
            if len(st) == 0:
                unpaired.append(letter)
            else:
                pairings.append((letter, min(st)))
                unpaired_b.remove(min(st))
        if len(unpaired_b) == 0:
            return None
        a = max(unpaired_b)
        s = 0
        while a + s + 1 in row_ip1:
            s += 1

        if a + s < row:
            return None
        new_row_ip1 = [let for let in row_ip1 if let != a]
        new_row_i = sorted([a + s, *row_i], reverse=True)
        ret_rc = RCGraph([*rc[: row - 1], tuple(new_row_i), tuple(new_row_ip1), *rc[row + 1 :]])
        if ret_rc.perm != rc.perm:
            return None
        
        ret = RootTableau.root_insert_rsk(ret_rc.perm_word, ret_rc.compatible_sequence)
        assert ret.edelman_greene_invariant == self.edelman_greene_invariant, f"{ret.edelman_greene_invariant=} != {self.edelman_greene_invariant=}"
        retmap = RootTableau.root_insert_rsk(self.reduced_word, self.compatible_sequence)

        # # mapper = {retmap.order_grid[ind]: ret.order_grid[ind] for ind in ret.iter_boxes}
        # box_list = []
        # for box in retmap.iter_boxes_row_word_order():
        #     root = retmap[box][0]
        #     # find root
        #     box2 = next(iter([box3 for box3 in self.iter_boxes if self[box3][0] == root]))
        #     box_list.append(box2)

        did = True
        while did:
            did = False
            for _box in ret.iter_outer_corners:
                if self._root_grid[_box] is not None:
                    ret = ret.up_jdt_slide(*_box, check=True)
                    did = True
        correct_word = _plactic_raising_operator(self.row_word, i)
        try_grid = copy.deepcopy(self._root_grid)
        # if ret.root_row_word != self.root_row_word:
        #     print("Contrast")
        #     print(f"{ret.root_row_word=}")
        #     print(f"{self.root_row_word=}")
        #     pretty_print(self)
        #     input()
        # relate order grid word grid
        for index, box in enumerate(self.iter_boxes_row_word_order()):
            try_grid[box] = (ret.perm.right_root_at(ret.order_grid[box], word=ret_rc.reduced_word), correct_word[index])
        retret = RootTableau(try_grid)
        return retret
        # for ind in self.iter_boxes:
        #     try_grid[ind] = (self.perm.right_root_at(self.order_grid[ind], word=ret.reduced_word), ret.compatible_sequence[self.order_grid[ind]])
        # assert ret == RootTableau(try_grid), "raising_operator construction mismatch"
        
        
        # for ind in np.ndindex(self._root_grid.shape):
        #     if self._root_grid[ind] is not None:

        #         # ind2 = np.where(retmap.order_grid == self.order_grid[ind])
        #         # print(f"Found that {self.order_grid[ind]=} is at")
        #         # print(ind2)
        #         # print("In")
        #         # print(f"{retmap.order_grid=}")
        #         # print(ret.order_grid[*ind2])
        #         try_grid[ind] = (ret.perm.right_root_at(self.order_grid[ind], word=ret.reduced_word), ret.compatible_sequence[self.order_grid[ind]])
        # tret = RootTableau(try_grid)
        # if any(self._root_grid[_box] is not None and ret._root_grid[_box] !=try_grid[_box] for _box in ret.iter_boxes):
        #     print("FYI ret")
        #     pretty_print(ret)
        #     print("tret")
        #     pretty_print(tret)
        #     print("self")
        #     pretty_print(self)
        #     input()
        
        
        
        #return ret
        # return ret

    def lowering_operator(self, row):
        # RF word is just the RC word backwards
        rc = self.rc_graph
        if row >= len(rc):
            return None
        row_i = [*rc[row - 1]]
        row_ip1 = [*rc[row]]

        # pair the letters
        pairings = []
        unpaired = []
        unpaired_b = [*row_ip1]

        for letter in row_i:
            st = [letter2 for letter2 in unpaired_b if letter2 > letter]
            if len(st) == 0:
                unpaired.append(letter)
            else:
                pairings.append((letter, min(st)))
                unpaired_b.remove(min(st))
        if len(unpaired) == 0:
            return None
        b = min(unpaired)
        t = min([j for j in range(b) if b - j - 1 not in row_i])

        if b - t < row + 1:
            return None
        new_row_i = [s for s in row_i if s != b]
        new_row_ip1 = sorted([b - t, *row_ip1], reverse=True)
        ret_rc = RCGraph([*rc[: row - 1], tuple(new_row_i), tuple(new_row_ip1), *rc[row + 1 :]])
        if ret_rc.perm != rc.perm:
            return None
        compatible_seq = ret_rc.compatible_sequence
        ret = RootTableau.root_insert_rsk(ret_rc.perm_word, compatible_seq)
        assert ret.edelman_greene_invariant == self.edelman_greene_invariant, f"{ret.edelman_greene_invariant=} != {self.edelman_greene_invariant=}"
        did = True
        while did:
            did = False
            for _box in ret.iter_outer_corners:
                if self._root_grid[_box] is not None:
                    ret = ret.up_jdt_slide(*_box, check=True)
                    did = True
        return ret

    @property
    def rc_graph(self):
        """
        Build an RCGraph from the tableau. If anything goes wrong during
        construction, raise a RuntimeError with a helpful debug payload
        (grid shape, a small snapshot of the grid, partial reduced/compatible
        sequences when available, and the original traceback).
        """
        reduced_word = None
        compatible_seq = None
        try:
            reduced_word, compatible_seq = _word_from_grid(self._root_grid, with_compatible_seq=True)
            rows = []
            if len(reduced_word) != len(compatible_seq):
                raise ValueError(f"length mismatch: reduced_word={len(reduced_word)} compatible_seq={len(compatible_seq)}")
            for i, a in enumerate(compatible_seq):
                if a > len(rows):
                    while len(rows) < a:
                        rows.append(())
                rows[a - 1] = (*rows[a - 1], reduced_word[i])
            return RCGraph(tuple(rows)).normalize()
        except Exception as exc:
            import traceback

            grid_shape = getattr(self._root_grid, "shape", None)
            grid_snapshot = repr(getattr(self, "_root_grid", None))

            tb = traceback.format_exc()
            raise RuntimeError(
                "Failed to build RCGraph from RootTableau.\n"
                f"grid.shape = {grid_shape}\n"
                f"grid_snapshot = {grid_snapshot}\n"
                f"reduced_word (partial) = {reduced_word}\n"
                f"compatible_seq (partial) = {compatible_seq}\n"
                f"traceback:\n{tb}"
            ) from exc

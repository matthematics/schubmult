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
    up_exists = (i - 1 >= 0 and (i - 1) < rows and j < cols and grid[i - 1, j] is not None)
    left_exists = (j - 1 >= 0 and (j - 1) < cols and i < rows and grid[i, j - 1] is not None)
    # consider hole on or beyond boundary as "outer"
    on_or_extends_boundary = (i >= rows) or (j >= cols) or (i == rows - 1) or (j == cols - 1)
    return (up_exists and left_exists) or ((up_exists or left_exists) and on_or_extends_boundary)


def _is_valid_inner_corner(grid: np.ndarray, i: int, j: int) -> bool:
    """
    Inner-corner predicate used by down_jdt_slide.
    Valid when the hole is inside the grid (must not extend grid) and there
    is a box below or to the right with the same boundary rules mirrored.
    """
    rows, cols = grid.shape
    if not (0 <= i < rows and 0 <= j < cols):
        return False
    if grid[i, j] is not None:
        return False
    down_exists = (i + 1 < rows and j < cols and grid[i + 1, j] is not None)
    right_exists = (j + 1 < cols and i < rows and grid[i, j + 1] is not None)
    on_boundary_up = (i == 0) or (j == 0)
    return (down_exists and right_exists) or ((down_exists or right_exists) and on_boundary_up)


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
                raise RuntimeError(f"_validate_grid: bad cell at {(i,j)}: {repr(cell)}")
            root, letter = cell[0], cell[1]
            try:
                int(letter)
            except Exception:
                raise RuntimeError(f"_validate_grid: non-int letter at {(i,j)}: {repr(letter)}")
            if not (isinstance(root, (tuple, list, int))):
                raise RuntimeError(f"_validate_grid: unexpected root at {(i,j)}: {repr(root)}")


def _snap_grid(grid: np.ndarray):
    """Return a compact, JSON-like snapshot of the grid for debug messages."""
    try:
        return [
            [None if grid[i, j] is None else (grid[i, j][0], int(grid[i, j][1])) for j in range(grid.shape[1])]
            for i in range(grid.shape[0])
        ]
    except Exception:
        return repr(grid)


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

    @property
    def perm(self):
        return Permutation.ref_product(*self.reduced_word)

    def rectify(self):
        cur = self
        while True:
            inner_corners = tuple(cur.iter_inner_corners())
            if not inner_corners:
                break
            cur = cur.down_jdt_slide(*next(iter(inner_corners)))
        return cur

    def up_jdt_slide(self, row, col):
        # capture preconditions & words
        try:
            _validate_grid(self._root_grid)
        except Exception as e:
            raise RuntimeError(f"up_jdt_slide: invalid input grid before slide: {e!r}")

        before_weight_tableau = self.weight_tableau
        before_reduced_word = self.reduced_word

        # allow holes that extend the grid by one; build an extended working grid
        grid = copy.deepcopy(self._root_grid)
        rows, cols = grid.shape
        if row >= rows or col >= cols:
            nr = max(rows, row + 1)
            nc = max(cols, col + 1)
            ext = np.empty((nr, nc), dtype=object)
            ext.fill(None)
            for r in range(rows):
                for c in range(cols):
                    ext[r, c] = grid[r, c]
            grid = ext
            rows, cols = grid.shape

        if not _is_valid_outer_corner(grid, row, col):
            raise ValueError("up_jdt_slide can only be performed from an outer-corner hole")

        # separate root-grid and letter-grid; letters move with the roots during jdt
        root_grid = np.empty(grid.shape, dtype=object)
        letter_grid = np.zeros(grid.shape, dtype=np.int32)
        for r0, c0 in np.ndindex(grid.shape):
            cell = grid[r0, c0]
            if cell is None:
                root_grid[r0, c0] = None
                letter_grid[r0, c0] = 0
            else:
                root_grid[r0, c0] = cell[0]
                letter_grid[r0, c0] = cell[1]

        r, c = row, col
        while True:
            above = (r - 1 >= 0 and root_grid[r - 1, c] is not None)
            left = (c - 1 >= 0 and root_grid[r, c - 1] is not None)
            if not above and not left:
                break
            if above and not left:
                # move from above into hole
                root_grid[r, c] = root_grid[r - 1, c]
                root_grid[r - 1, c] = None
                letter_grid[r, c] = letter_grid[r - 1, c]
                letter_grid[r - 1, c] = 0
                r -= 1
                continue
            if left and not above:
                # move from left into hole
                root_grid[r, c] = root_grid[r, c - 1]
                root_grid[r, c - 1] = None
                letter_grid[r, c] = letter_grid[r, c - 1]
                letter_grid[r, c - 1] = 0
                c -= 1
                continue
            # both available: choose by letter_grid (recording letter) first
            val_above = letter_grid[r - 1, c]
            val_left = letter_grid[r, c - 1]
            
            if val_above >= val_left:
                root_grid[r, c] = root_grid[r - 1, c]
                root_grid[r - 1, c] = None
                letter_grid[r, c] = val_above
                letter_grid[r - 1, c] = 0
                r -= 1
            else:
                root_grid[r, c] = root_grid[r, c - 1]
                root_grid[r, c - 1] = None
                letter_grid[r, c] = val_left
                letter_grid[r, c - 1] = 0
                c -= 1
        root_grid[row, col] = None
        letter_grid[row, col] = 0
        # construct combined grid
        new_grid = np.empty(root_grid.shape, dtype=object)
        for r0, c0 in np.ndindex(root_grid.shape):
            if root_grid[r0, c0] is None:
                new_grid[r0, c0] = None
            else:
                new_grid[r0, c0] = (root_grid[r0, c0], int(letter_grid[r0, c0]))

        try:
            _validate_grid(new_grid)
        except Exception as e:
            raise RuntimeError(
                f"up_jdt_slide produced invalid grid after slide: {e!r}\nfinal_hole={(r,c)}\ngrid_snapshot={_snap_grid(new_grid)}"
            )

        after = RootTableau(new_grid)
        # sanity-check invariants: plactic (weight) tableau should be unchanged
        if self.rc_graph != after.rc_graph:
            raise ValueError(f"Does not preserve bunkbaby {self.rc_graph=} {after.rc_graph=}")
        
        return after

    def down_jdt_slide(self, row, col):
        """
        Perform a downward/rightward jeu-de-taquin slide starting from the given
        (row, col) hole (0-indexed). Boxes from below or to the right are moved
        into the hole, preferring the smaller recording letter when both exist.
        Returns a new RootTableau (does not mutate self).
        """
        try:
            _validate_grid(self._root_grid)
        except Exception as e:
            raise RuntimeError(f"down_jdt_slide: invalid input grid before slide: {e!r}")

        if not _is_valid_inner_corner(self._root_grid, row, col):
            raise ValueError("down_jdt_slide can only be performed from an inner-corner hole")

        # separate root and letter grids; letters move with the roots during jdt
        root_grid = np.empty(self._root_grid.shape, dtype=object)
        letter_grid = np.empty(self._root_grid.shape, dtype=object)
        for r0, c0 in np.ndindex(self._root_grid.shape):
            cell = self._root_grid[r0, c0]
            if cell is None:
                root_grid[r0, c0] = None
                letter_grid[r0, c0] = None
            else:
                root_grid[r0, c0] = cell[0]
                letter_grid[r0, c0] = cell[1]

        before_weight_tableau = self.weight_tableau
        before_reduced_word = self.reduced_word

        new_root_grid = copy.deepcopy(root_grid)
        new_letter_grid = copy.deepcopy(letter_grid)
        rows, cols = new_root_grid.shape
        r, c = row, col

        while True:
            down_exists = (r + 1 < rows and new_root_grid[r + 1, c] is not None)
            right_exists = (c + 1 < cols and new_root_grid[r, c + 1] is not None)
            if not down_exists and not right_exists:
                break
            if down_exists and not right_exists:
                new_root_grid[r, c] = new_root_grid[r + 1, c]
                new_root_grid[r + 1, c] = None
                new_letter_grid[r, c] = new_letter_grid[r + 1, c]
                new_letter_grid[r + 1, c] = None
                r += 1
                continue
            if right_exists and not down_exists:
                new_root_grid[r, c] = new_root_grid[r, c + 1]
                new_root_grid[r, c + 1] = None
                new_letter_grid[r, c] = new_letter_grid[r, c + 1]
                new_letter_grid[r, c + 1] = None
                c += 1
                continue
            # both present: choose by anchored recording letter
            down_letter = new_letter_grid[r + 1, c]
            right_letter = new_letter_grid[r, c + 1]
            if down_letter <= right_letter:
                new_root_grid[r, c] = new_root_grid[r + 1, c]
                new_root_grid[r + 1, c] = None
                new_letter_grid[r, c] = new_letter_grid[r + 1, c]
                new_letter_grid[r + 1, c] = None
                r += 1
            else:
                new_root_grid[r, c] = new_root_grid[r, c + 1]
                new_root_grid[r, c + 1] = None
                new_letter_grid[r, c] = new_letter_grid[r, c + 1]
                new_letter_grid[r, c + 1] = None
                c += 1

        # assemble final grid by attaching moved letters back to their moved roots
        new_grid = np.empty(new_root_grid.shape, dtype=object)
        for r0, c0 in np.ndindex(new_root_grid.shape):
            if new_root_grid[r0, c0] is None:
                new_grid[r0, c0] = None
            else:
                new_grid[r0, c0] = (new_root_grid[r0, c0], int(new_letter_grid[r0, c0]) if new_letter_grid[r0, c0] is not None else None)

        try:
            _validate_grid(new_grid)
        except Exception as e:
            raise RuntimeError(
                f"down_jdt_slide produced invalid grid after slide: {e!r}\nfinal_hole={(r,c)}\ngrid_snapshot={_snap_grid(new_grid)}"
            )

        after = RootTableau(new_grid)
        # sanity-check plactic invariant
        if before_weight_tableau != after.weight_tableau or before_reduced_word != after.reduced_word:
            raise RuntimeError(
                "down_jdt_slide violated invariants\n"
                f"before_weight_tableau={before_weight_tableau} after_weight_tableau={after.weight_tableau}\n"
                f"before_reduced_word={before_reduced_word} after_reduced_word={after.reduced_word}\n"
                f"before_grid={_snap_grid(self._root_grid)}\nafter_grid={_snap_grid(new_grid)}"
            )

        return after

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

    def iter_boxes(self):
        for i in range(self.rows):
            for j in range(self.cols):
                cell = self._root_grid[i, j]
                if cell is not None:
                    yield (i, j)

    def iter_outer_corners(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if _is_valid_outer_corner(self._root_grid, i, j):
                    yield (i, j)

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
        rows = []
        for row in self._root_grid:
            rows.append([])
            for cell in row:
                if cell is None:
                    rows[-1].append(0)
                else:
                    rows[-1].append(cell[1])

        return Plactic(rows)

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

    def raising_operator(self, i):
        print("Should work for non-jdt cases")
        rc = self.rc_graph.raising_operator(i)
        if rc is None:
            return None
        return RootTableau.from_rc_graph(rc)
        # """Crystal raising operator e_i on the root tableau"""
        # word = [*self.row_word]
        # opening_stack = []
        # closing_stack = []
        # for index in range(len(word)):
        #     if word[index] == i + 1:
        #         opening_stack.append(index)
        #     elif word[index] == i:
        #         if len(opening_stack) > 0:
        #             opening_stack.pop()
        #         else:
        #             closing_stack.append(index)
        # if len(opening_stack) == 0:
        #     return None
        # index_to_change = opening_stack[0]
        # word[index_to_change] = i
        # new_grid = copy.deepcopy(self._root_grid)
        # the_index = 0
        # for r in range(self.rows - 1, -1, -1):
        #     for c in range(self.cols):
        #         cell = self._root_grid[r, c]
        #         if cell is not None:
        #             if the_index == index_to_change:
        #                 root_cell, _ = cell
        #                 new_grid[r, c] = (root_cell, i)
        #                 print("This doesn't work")
        #                 return RootTableau(new_grid)
        #             the_index += 1
        # return None

        # # RF word is just the RC word backwards
        # if row >= len(self):
        #     return None
        # row_i = [*self[row - 1]]
        # row_ip1 = [*self[row]]

        # # pair the letters
        # pairings = []
        # unpaired = []
        # unpaired_b = [*row_ip1]

        # for letter in row_i:
        #     st = [letter2 for letter2 in unpaired_b if letter2 > letter]
        #     if len(st) == 0:
        #         unpaired.append(letter)
        #     else:
        #         pairings.append((letter, min(st)))
        #         unpaired_b.remove(min(st))
        # if len(unpaired_b) == 0:
        #     return None
        # a = max(unpaired_b)
        # s = 0
        # while a + s + 1 in row_ip1:
        #     s += 1

        # if a + s < row:
        #     return None
        # new_row_ip1 = [let for let in row_ip1 if let != a]
        # new_row_i = sorted([a + s, *row_i], reverse=True)
        # ret_rc = type(self)([*self[: row - 1], tuple(new_row_i), tuple(new_row_ip1), *self[row + 1 :]])
        # if ret_rc.perm != self.perm:
        #     return None
        # return ret_rc

    def lowering_operator(self, index):
        print("Should work for non-jdt cases")
        rc = self.rc_graph.lowering_operator(index)
        if rc is None:
            return None
        return RootTableau.from_rc_graph(rc)

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

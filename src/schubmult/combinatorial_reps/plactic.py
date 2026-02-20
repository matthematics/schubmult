# filepath: /home/matthematics/schubmult/src/schubmult.combinatorial_reps.plactic.py
import copy

import numpy as np

from schubmult.combinatorial_reps.permutation import Permutation  # noqa: F401

from ..utils._grid_print import GridPrint
from .crystal_graph import CrystalGraph


def _is_valid_outer_corner(grid, i: int, j: int) -> bool:
    """Outer-corner predicate: empty position where BOTH up and left have filled cells or are at boundary.
    Formal: position where we can place a hole, with both (i-1,j) and (i,j-1) in the shape."""
    rows, cols = grid.shape

    # Position (0,0) is never an outer corner (no up or left neighbors)
    if i == 0 and j == 0:
        return False

    # Position must be empty if within grid
    if 0 <= i < rows and 0 <= j < cols:
        if grid[i, j] is not None:
            return False

    # BOTH up and left must have content (be in the shape)
    # If we're at i=0 or j=0, that direction is at boundary (counts as having content)
    up_valid = False
    left_valid = False

    if i == 0:
        # At top boundary, up is valid
        up_valid = True
    elif i - 1 < rows and j < cols:
        up_valid = grid[i - 1, j] is not None

    if j == 0:
        # At left boundary, left is valid
        left_valid = True
    elif i < rows and j - 1 < cols:
        left_valid = grid[i, j - 1] is not None

    return up_valid and left_valid
    # consider hole on or beyond boundary as "outer


def _is_valid_inner_corner(grid, i: int, j: int, inner_shape=None) -> bool:
    """Inner-corner predicate: a hole (position in inner_shape) where BOTH down and right are NOT holes.

    An inner corner is a cell that IS in the inner shape (a hole) where you cannot slide it
    down-right any further because both the cell below and the cell to the right are either
    filled (have values) or outside the inner shape (but may be empty cells in outer shape).
    """
    rows, cols = grid.shape

    # Must be within grid bounds
    if i < 0 or j < 0 or i >= rows or j >= cols:
        return False

    # If no inner_shape, use the old logic: position must be empty (None or 0)
    if inner_shape is None:
        if grid[i, j] is not None and grid[i, j] != 0:
            return False

        # Check down neighbor: must be outside grid OR filled
        down_blocked = False
        if i + 1 >= rows:
            down_blocked = True
        elif grid[i + 1, j] is not None and grid[i + 1, j] != 0:
            down_blocked = True

        # Check right neighbor: must be outside grid OR filled
        right_blocked = False
        if j + 1 >= cols:
            right_blocked = True
        elif grid[i, j + 1] is not None and grid[i, j + 1] != 0:
            right_blocked = True

        return down_blocked and right_blocked

    # With inner_shape: position MUST be a hole (in inner_shape) AND actually empty
    if not (i < len(inner_shape) and j < inner_shape[i]):
        return False

    # Position must actually be empty (not filled in)
    if grid[i, j] is not None and grid[i, j] != 0:
        return False

    # Check down neighbor: must NOT be a hole (either filled, empty-in-outer-shape, or outside grid)
    down_not_hole = False
    down_has_value = False
    if i + 1 >= rows:
        # Beyond grid boundary - not a hole
        down_not_hole = True
    elif i + 1 >= len(inner_shape) or j >= inner_shape[i + 1]:
        # Outside inner shape (either filled or empty cell in outer shape) - not a hole
        down_not_hole = True
        if grid[i + 1, j] is not None and grid[i + 1, j] != 0:
            down_has_value = True
    elif grid[i + 1, j] is not None and grid[i + 1, j] != 0:
        # Inside inner_shape but has a value - not a hole anymore
        down_not_hole = True
        down_has_value = True

    # Check right neighbor: must NOT be a hole (either filled, empty-in-outer-shape, or outside grid)
    right_not_hole = False
    right_has_value = False
    if j + 1 >= cols:
        # Beyond grid boundary - not a hole
        right_not_hole = True
    elif i >= len(inner_shape) or j + 1 >= inner_shape[i]:
        # Outside inner shape (either filled or empty cell in outer shape) - not a hole
        right_not_hole = True
        if grid[i, j + 1] is not None and grid[i, j + 1] != 0:
            right_has_value = True
    elif grid[i, j + 1] is not None and grid[i, j + 1] != 0:
        # Inside inner_shape but has a value - not a hole anymore
        right_not_hole = True
        right_has_value = True

    # Inner corner only if BOTH directions are NOT holes AND at least one has a value
    return down_not_hole and right_not_hole and (down_has_value or right_has_value)


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


class Plactic(GridPrint, CrystalGraph):
    def __iter__(self):
        for i in range(self._grid.shape[0]):
            for j in range(self._grid.shape[1]):
                if self._grid[i, j] is not None:
                    yield self._grid[i, j]

    def __le__(self, other):
        if not isinstance(other, Plactic):
            return NotImplemented
        return self.row_word <= other.row_word

    def __lt__(self, other):
        if not isinstance(other, Plactic):
            return NotImplemented
        return self.row_word < other.row_word

    # in order of row row
    @property
    def iter_boxes(self):
        # Iterate over non-border cells only
        for i in range(self.rows - 1, -1, -1):
            for j in range(self.cols):
                if self._grid[i, j] is not None:
                    yield (i, j)

    def up_jdt_slide(self, row, col):
        """
        Perform an upward jeu de taquin slide starting from the given (row, col)
        position (0-indexed). Returns a new Plactic tableau.
        """
        new_grid = copy.deepcopy(self._grid)
        i, j = row, col
        # Resize grid if needed
        if self.rows <= i or self.cols <= j:
            new_grid.resize((max(self.rows, i + 1), max(self.cols, j + 1)), refcheck=False)
        if self[i, j] != 0 and self[i, j] is not None:
            raise ValueError(f"up_jdt_slide starting position must be empty, got {self[i, j]=}")
        if self[i - 1, j] in (0, None) and self[i, j - 1] in (0, None):
            raise ValueError("up_jdt_slide starting position has no valid moves")

        while True:
            up = new_grid[i - 1, j] if i - 1 >= 0 else None
            left = new_grid[i, j - 1] if j - 1 >= 0 else None

            if up in (None, 0) and left in (None, 0):
                # No more moves possible
                break
            if up in (None, 0):
                # Can only move left
                new_grid[i, j] = left
                j -= 1
            elif left in (None, 0):
                # Can only move up
                new_grid[i, j] = up
                i -= 1
            else:
                # Both moves possible; choose the larger entry
                if up >= left:
                    new_grid[i, j] = up
                    i -= 1
                else:
                    new_grid[i, j] = left
                    j -= 1

        new_grid[i, j] = 0
        new_inner_shape = [*self._inner_shape]
        new_inner_shape[i] += 1

        return Plactic._from_grid(new_grid, tuple(new_inner_shape))

    def down_jdt_slide(self, row, col):
        """
        Perform a jeu de taquin slide starting from the given (row, col)
        position (0-indexed). Returns a new Plactic tableau.
        """
        new_grid = copy.deepcopy(self._grid)
        i, j = row, col
        if self[i, j] != 0 and self[i, j] is not None:
            raise ValueError(f"down_jdt_slide starting position must be empty, got {self[i, j]=}")
        if self[i + 1, j] in (0, None) and self[i, j + 1] in (0, None):
            raise ValueError("down_jdt_slide starting position has no valid moves")
        while True:
            down = new_grid[i + 1, j] if i + 1 < self.rows else None
            right = new_grid[i, j + 1] if j + 1 < self.cols else None

            if down in (None, 0) and right in (None, 0):
                # No more moves possible
                break
            if down in (None, 0):
                # Can only move right
                new_grid[i, j] = right
                j += 1
            elif right in (None, 0):
                # Can only move down
                new_grid[i, j] = down
                i += 1
            else:
                # Both moves possible; choose the smaller entry
                if down <= right:
                    new_grid[i, j] = down
                    i += 1
                else:
                    new_grid[i, j] = right
                    j += 1

        # Mark the final position as None
        new_grid[i, j] = None

        new_inner_shape = [*self._inner_shape] # if we are doing down_jdt, inner shape had better not be None
        new_inner_shape[row] -= 1
        return Plactic._from_grid(new_grid, tuple(new_inner_shape))

    @property
    def iter_outer_corners(self):
        # Check border positions for outer corners
        # Check positions at (i, self.cols) for all rows and (self.rows, j) for all cols
        for i in range(self._grid.shape[0]):
            for j in range(self._grid.shape[1]):
                if _is_valid_outer_corner(self._grid, i, j):
                    yield (i, j)

    @property
    def iter_inner_corners(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if _is_valid_inner_corner(self._grid, i, j, self._inner_shape):
                    yield (i, j)

    def __init__(self, word=(), inner_shape=None):
        # Convert word (tuple of tuples) to np.ndarray
        # Grid always has +1 row and +1 col as empty border (ensures outer corners exist)
        if isinstance(word, np.ndarray):
            # Assume input grid already has border if it's an ndarray
            self._grid = word
            self._inner_shape = inner_shape
        elif len(word) == 0:
            # Empty tableau: just 1x1 border
            self._grid = np.full((1, 1), None, dtype=object)
            self._inner_shape = None
        else:
            # Find max row length
            max_cols = max((len(r) for r in word), default=0)
            num_rows = len(word)
            # Add +1 to both dimensions for border
            self._grid = np.full((num_rows + 1, max_cols + 1), None, dtype=object)
            for i, row in enumerate(word):
                for j, val in enumerate(row):
                    if val != 0:  # Don't store explicit zeros, use None for empty
                        self._grid[i, j] = int(val)
            # Track inner shape: infer from 0 values if not provided
            if inner_shape is None:
                inner_shape_list = []
                for i, row in enumerate(word):
                    # Count leading zeros/Nones in this row
                    count = 0
                    for val in row:
                        if val == 0 or val is None:
                            count += 1
                        else:
                            break
                    inner_shape_list.append(count)
                # Only set if there are actual holes
                if any(c > 0 for c in inner_shape_list):
                    self._inner_shape = tuple(inner_shape_list)
                else:
                    self._inner_shape = None
            else:
                self._inner_shape = tuple(int(x) for x in inner_shape)

    @classmethod
    def _from_grid(cls, grid, inner_shape=None):
        """Create a Plactic directly from an np.ndarray grid.
        Assumes grid already includes the border."""
        obj = cls.__new__(cls)
        obj._grid = grid
        obj._inner_shape = inner_shape
        return obj

    def shiftup(self, k):
        """Return a new Plactic with all entries increased by k."""
        new_grid = copy.deepcopy(self._grid)
        # Only process non-border cells
        for i in range(self.rows):
            for j in range(self.cols):
                if new_grid[i, j] is not None and new_grid[i, j] != 0:
                    new_grid[i, j] = int(new_grid[i, j]) + k
        return Plactic._from_grid(new_grid, self._inner_shape)

    @classmethod
    def all_ss_tableaux(cls, shape, max_entry, inner_shape=None):
        """Generate all semistandard tableaux of given shape (or skew shape) with entries <= max_entry.

        Args:
            shape: Sequence of row lengths (outer shape)
            max_entry: Maximum entry value
            inner_shape: Optional sequence of left offsets per row (for skew shapes).
                        If provided, positions [row][0:inner_shape[row]] are marked as 0.

        Returns:
            Set of Plactic instances representing all valid semistandard tableaux
        """
        # Normalize shapes
        outer = tuple(int(x) for x in shape)
        r = len(outer)
        if inner_shape is None:
            inner = (0,) * r
        else:
            inner = tuple(int(x) for x in inner_shape) + (0,) * (r - len(inner_shape))

        tableaux = set()

        # Initialize tableau with 0s for inner cells
        current_tableau = []
        for i in range(r):
            row = [0] * outer[i]
            current_tableau.append(row)

        # Build list of cells to fill (excluding inner cells)
        cells = []
        for i in range(r):
            for j in range(inner[i], outer[i]):
                cells.append((i, j))

        def _recurse(cell_idx):
            if cell_idx == len(cells):
                # Base case: A complete tableau is found
                tabby = cls([tuple(r) for r in current_tableau], inner_shape=inner if any(inner) else None)
                tabby._grid[tabby._grid == 0] = None  # Convert explicit zeros to None
                tableaux.add(tabby)
                return

            row, col = cells[cell_idx]

            # Determine the minimum possible value for the current cell
            # Semistandard: rows weakly increasing, columns strictly increasing
            min_val = 1

            # Check left neighbor (same row) - must be >= (weakly increasing rows)
            if col > 0:
                left_val = current_tableau[row][col - 1]
                if left_val != 0:  # Skip inner cells
                    min_val = max(min_val, left_val)

            # Check top neighbor (same column, previous row) - must be > (strictly increasing columns)
            if row > 0 and col < len(current_tableau[row - 1]):
                top_val = current_tableau[row - 1][col]
                if top_val != 0:  # Skip inner cells
                    min_val = max(min_val, top_val + 1)

            for val in range(min_val, max_entry + 1):
                current_tableau[row][col] = val
                _recurse(cell_idx + 1)
                # Backtrack (set back to 0 for clarity)
                current_tableau[row][col] = 0

        _recurse(0)
        return tableaux

    @property
    def row_word(self):
        """Return the row-reading word as a flat tuple."""
        word = []
        # Read from bottom to top, excluding border row
        for i in range(self.rows - 1, -1, -1):
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    word.append(val)
        return tuple(word)

    @property
    def column_word(self):
        wrd = []
        # Exclude border column
        for j in range(self.cols):
            for i in range(self.rows - 1, -1, -1):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    wrd.append(val)
        return tuple(wrd)

    def transpose(self):
        """Return the transpose of this Plactic tableau."""
        if self.rows == 0 or self.cols == 0:
            return self.__class__(())
        # Create transposed grid with border (cols+1 by rows+1)
        new_grid = np.full((self.cols + 1, self.rows + 1), None, dtype=object)
        for i in range(self.rows):
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None:
                    new_grid[j, i] = val

        # Transpose inner_shape: convert row offsets to column offsets
        inner_shape = self._inner_shape
        transposed_inner = None
        if inner_shape is not None:
            # Build transposed inner shape
            transposed_inner = [0] * self.cols
            for i in range(len(inner_shape)):
                for j in range(inner_shape[i]):
                    if j < len(transposed_inner):
                        transposed_inner[j] += 1
            transposed_inner = tuple(transposed_inner) if any(x > 0 for x in transposed_inner) else None

        return self.__class__._from_grid(new_grid, transposed_inner)

    @property
    def rows(self):
        # Exclude the border row
        return max(0, self._grid.shape[0] - 1)

    @property
    def cols(self):
        # Exclude the border column
        return max(0, self._grid.shape[1] - 1)

    def __getitem__(self, key):
        # FLIPPED FOR PRINTING
        if isinstance(key, int):
            # Return row as tuple for compatibility
            return tuple(self._grid[key, j] for j in range(self._grid.shape[1]))
        if isinstance(key, tuple):
            i, j = key
            if i >= self._grid.shape[0] or j >= self._grid.shape[1] or i < 0 or j < 0:
                return None
            return self._grid[i, j]
        is_slice = isinstance(key, slice)

        if is_slice:
            return tuple(tuple(self._grid[n, j] for j in range(self._grid.shape[1])) for n in range(self.rows)[key])

        raise ValueError(f"Bad indexing {key=}")

    @staticmethod
    def _remap_word(word, fn):
        """Apply fn to every entry of a tableau-like tuple-of-rows."""
        return tuple(tuple(fn(x) for x in row) for row in tuple(word))

    def invert(self):
        """
        Return a Plactic whose entries are remapped so that standard
        (increasing) insertion order applies. If reverse_semistandard is True
        we negate entries (so larger original becomes smaller).
        """
        new_grid = copy.deepcopy(self._grid)
        # Only invert non-border cells
        for i in range(self.rows):
            for j in range(self.cols):
                if new_grid[i, j] is not None and new_grid[i, j] != 0:
                    new_grid[i, j] = -int(new_grid[i, j])
        return Plactic._from_grid(new_grid, self._inner_shape)

    def __mul__(self, other):
        """
        Plactic product: insert entries of `other` in row-reading order
        (top-to-bottom, left-to-right) into a copy of self.
        """
        if not isinstance(other, Plactic):
            return NotImplemented
        word = [*self.row_word, *other.row_word]
        pl = Plactic()
        return pl.rs_insert(*word)

    @property
    def shape(self):
        shape_list = []
        # Only count non-border rows
        for i in range(self.rows):
            count = sum(1 for j in range(self.cols) if self._grid[i, j] is not None and self._grid[i, j] != 0)
            if count > 0:
                shape_list.append(count)
        return tuple(shape_list)

    @property
    def skew_shape(self):
        """Return the skew shape as a tuple of (row_length, left_offset) pairs."""
        outer_shape = []
        inner_shape = []
        for i in range(self.rows):
            row_length = 0
            left_offset = 0
            for j in range(self.cols):
                if self._grid[i, j] is not None and self._grid[i, j] != 0:
                    row_length += 1
                elif row_length == 0:
                    left_offset += 1
            if row_length > 0:
                outer_shape.append(left_offset + row_length)
                #shape_list.append((row_length, left_offset))
                inner_shape.append(left_offset)
        return tuple(outer_shape), tuple(inner_shape)

    @classmethod
    def from_word(cls, word):
        # accept any iterable-of-rows and normalize to tuple-of-tuples
        return cls().rs_insert(*word)

    def __hash__(self):
        # Hash based on the actual tableau content (excluding border)
        content = []
        for i in range(self.rows):
            row_vals = []
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    row_vals.append(val)
            if row_vals:
                content.append(tuple(row_vals))
        return hash(tuple(content))

    @staticmethod
    def _rs_insert(word, letter, i=0):
        """
        Row-insertion: insert `letter` into `word` (tuple-of-tuples) starting at row i.
        Returns a new tuple-of-tuples representing the tableau after insertion.

        Algorithm:
        - If row i doesn't exist, append a new row with the letter.
        - Otherwise find the first entry in row i that is > letter, replace it (bump)
          and recursively insert the bumped value into row i+1.
        - If no entry > letter, append letter to the end of row i.
        """
        # normalize input
        word = tuple(tuple(r) for r in word)

        if i >= len(word):
            # append new row containing the letter
            return (*word, (int(letter),))

        row = list(word[i])
        x = int(letter)

        # find first entry strictly greater than x to bump
        for k, val in enumerate(row):
            if val is not None and val != 0 and val > x:
                bumped = val
                row[k] = x
                new_word = (*word[:i], tuple(row), *word[i + 1 :])
                # recursively insert bumped into next row
                return Plactic._rs_insert(new_word, bumped, i=i + 1)

        # nothing to bump -> append to end of this row
        row.append(x)
        return (*word[:i], tuple(row), *word[i + 1 :])

    def rs_insert(self, *letters):
        """Insert one or more letters in sequence (row-insertion) and return a new Plactic."""
        # accept either rs_insert(a,b,...) or rs_insert([a,b,...])
        if len(letters) == 1 and isinstance(letters[0], (list, tuple)):
            seq = list(letters[0])
        else:
            seq = list(letters)

        # Convert current grid to tuple-of-tuples format for insertion (excluding border)
        word_tuples = []
        for i in range(self.rows):
            row_vals = []
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    row_vals.append(val)
            if row_vals:  # Only add non-empty rows
                word_tuples.append(tuple(row_vals))
        new_word = tuple(word_tuples)

        for letter in seq:
            new_word = Plactic._rs_insert(new_word, int(letter))
        return Plactic(new_word)

    def __eq__(self, other):
        if isinstance(other, Plactic):
            newfangled = self.__class__.from_word(self.row_word)
            oldfangled = self.__class__.from_word(other.row_word)
            if newfangled.row_word != oldfangled.row_word:
                return False
            return True
        return False

    # -- CrystalGraph API wrappers (delegate to the NilPlactic <-> RCGraph machinery) --
    # These adapt the RCGraph crystal operators/queries to Plactic tableaux by
    # converting the tableau to its highest-weight RC graph (via NilPlactic.hw_rc),
    # delegating to the RCGraph implementation, then mapping results back to Plactic.

    def raising_operator(self, i):
        """Crystal raising operator e_i on the Plactic tableau (delegates to RCGraph)."""
        word = [*self.row_word]
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
        result = Plactic().rs_insert(*word)
        assert result.shape == self.shape, f"{result.shape=} {self.shape=}"
        return result

    def lowering_operator(self, i):
        """Crystal lowering operator f_i on the Plactic tableau (delegates to RCGraph)."""
        word = [*self.row_word]
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
        if len(closing_stack) == 0:
            return None
        index_to_change = closing_stack[-1]
        word[index_to_change] = i + 1
        result = Plactic().rs_insert(*word)
        assert result.shape == self.shape, f"{result.shape=} {self.shape=} {result=} {self=}"
        return result

    # def _bracket_positions(self, i: int):
    #     """
    #     Bracketing scan for indices i (as '+') and i+1 (as '-') on the tableau's
    #     row-reading word. Returns (plus_stack, minus_unpaired, word_list) where:
    #       - plus_stack is a list of indices of unmatched i (these are the '+' positions)
    #       - minus_unpaired is a list of indices of unmatched i+1 (these are the '-' positions)
    #       - word_list is a mutable list copy of the row-reading word used for modification.
    #     The scan follows the standard left->right cancellation rule.
    #     """
    #     word = list(self.row_word)
    #     plus_stack = []
    #     minus_unpaired = []
    #     for idx, a in enumerate(word):
    #         if a == i:
    #             plus_stack.append(idx)
    #         elif a == i + 1:
    #             if plus_stack:
    #                 # cancel with the most recent unmatched '+'
    #                 plus_stack.pop()
    #             else:
    #                 # unmatched '-'
    #                 minus_unpaired.append(idx)
    #     return plus_stack, minus_unpaired, word

    # def raising_operator(self, i):
    #     """
    #     Crystal raising operator e_i: change the leftmost unmatched (i+1) -> i,
    #     if any; rebuild the Plactic from the modified row-reading word and return it.
    #     Returns None if epsilon_i == 0.
    #     """
    #     plus, minus, word = self._bracket_positions(i)
    #     if len(minus) == 0:
    #         return None
    #     index_to_change = minus[0]  # leftmost unmatched i+1
    #     word[index_to_change] = i
    #     result = Plactic().rs_insert(*word)
    #     # shape-preserving check (crystal operators preserve shape)
    #     assert result.shape == self.shape, f"shape changed by e_{i}: {self.shape} -> {result.shape}"
    #     return result

    # def lowering_operator(self, i):
    #     """
    #     Crystal lowering operator f_i: change the rightmost unmatched i -> i+1,
    #     if any; rebuild the Plactic from the modified row-reading word and return it.
    #     Returns None if phi_i == 0.
    #     """
    #     plus, minus, word = self._bracket_positions(i)
    #     if len(plus) == 0:
    #         return None
    #     index_to_change = plus[-1]  # rightmost unmatched i
    #     word[index_to_change] = i + 1
    #     result = Plactic().rs_insert(*word)
    #     assert result.shape == self.shape, f"shape changed by f_{i}: {self.shape} -> {result.shape}"
    #     return result

    @property
    def crystal_weight(self):
        """Return the crystal weight of this tableau (delegated to RCGraph)."""
        wt = []
        for i in range(1, max(self.row_word, default=0) + 1):
            wt.append(sum(1 for a in self.row_word if a == i))
        return tuple(wt)

    def crystal_length(self):
        """Return the length/number of rows used for the crystal"""
        return 100000000000000

    @classmethod
    def yamanouchi(cls, shape):
        """
        Return the Yamanouchi (highest-weight) tableau of the given shape.
        """
        if not shape:
            return cls(())
        num_rows = len(shape)
        max_cols = max(shape) if shape else 0
        # Create grid with border
        grid = np.full((num_rows + 1, max_cols + 1), None, dtype=object)
        for i in range(num_rows):
            for j in range(shape[i]):
                grid[i, j] = i + 1
        return cls._from_grid(grid)

    @property
    def is_increasing(self):
        """Check if the tableau is strictly increasing in rows and columns."""
        for i in range(self.rows):
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is None or val == 0:
                    continue
                # Check row condition: increasing left to right
                if j > 0:
                    left_val = self._grid[i, j - 1]
                    if left_val is not None and left_val != 0 and val <= left_val:
                        return False
                # Check column condition: increasing top to bottom
                if i > 0:
                    above_val = self._grid[i - 1, j]
                    if above_val is not None and above_val != 0 and val <= above_val:
                        return False
        return True

    def rectify(self):
        # if self.rows == 0:
        #     return self
        # if self.cols == 0:
        #     return self

        # # Check if first row's first element is 0 or None
        # first_val = self._grid[0, 0]
        # if first_val is not None and first_val != 0:
        #     return self

        # # Find first non-zero/non-None element in first row
        # index = 0
        # while index < self.cols and (self._grid[0, index] in (0, None)):
        #     index += 1

        # if index == self.cols:
        #     # Entire first row is empty
        #     if self.rows == 1:
        #         return self
        #     # Check if second row starts with 0/None
        #     if self.rows > 1 and (self._grid[1, 0] in (0, None)):
        #         # Recursively rectify remaining rows
        #         remaining_grid = self._grid[1:, :]
        #         remaining = self.__class__._from_grid(remaining_grid).rectify()
        #         # Prepend the empty first row
        #         result_grid = np.vstack([self._grid[0:1, :], remaining._grid]) if remaining.rows > 0 else self._grid[0:1, :]
        #         return self.__class__._from_grid(result_grid)
        #     return self.down_jdt_slide(0, 0).rectify()

        # index -= 1
        # return self.down_jdt_slide(0, index).rectify()
        ret = self.__class__._from_grid(self._grid.copy(), self._inner_shape)
        while True:
            try:
                corner = next(ret.iter_inner_corners)
                ret = ret.down_jdt_slide(*corner)
            except StopIteration:
                break
        return ret

    @classmethod
    def superstandard(cls, shape):
        if shape is None:
            return None
        if not shape:
            return cls(())
        num_rows = len(shape)
        max_cols = max(shape) if shape else 0
        # Create grid with border
        grid = np.full((num_rows + 1, max_cols + 1), None, dtype=object)
        index = 1
        for i in range(num_rows):
            for j in range(shape[i]):
                grid[i, j] = index
                index += 1
        return cls._from_grid(grid)

    @property
    def is_semistandard(self):
        # Only check non-border cells
        for i in range(self.rows):
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is None or val == 0:
                    continue
                # Check row condition: increasing left to right
                if j > 0:
                    left_val = self._grid[i, j - 1]
                    if left_val is not None and left_val != 0 and val < left_val:
                        return False
                # Check column condition: strictly increasing top to bottom
                if i > 0:
                    above_val = self._grid[i - 1, j]
                    if above_val is not None and above_val != 0 and val <= above_val:
                        return False
        return True

    def reverse_rsk(self, recording_tableau):
        """
        Inverse RSK (row-insertion) for the pair (P,Q) where `self` is P and
        `recording_tableau` is the standard recording tableau Q of the same shape.

        Returns the original word as a list of integers (in insertion order).
        """
        # Create mutable copies by converting grids to list of lists (excluding border)
        P = []
        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None and val != 0:
                    row.append(val)
            if row:
                P.append(row)

        Q = []
        for i in range(recording_tableau.rows):
            row = []
            for j in range(recording_tableau.cols):
                val = recording_tableau._grid[i, j]
                if val is not None and val != 0:
                    row.append(val)
            if row:
                Q.append(row)

        total = sum(len(r) for r in Q)
        word_rev = []

        for t in range(total, 0, -1):
            # find position (i,j) of t in Q
            pos_i = pos_j = None
            for i, row in enumerate(Q):
                for j, val in enumerate(row):
                    if val == t:
                        pos_i, pos_j = i, j
                        break
                if pos_i is not None:
                    break
            if pos_i is None:
                raise ValueError(f"recording tableau missing entry {t}")

            i, j = pos_i, pos_j

            # remove the box from Q
            Q[i].pop(j)
            if len(Q[i]) == 0:
                Q.pop(i)

            # remove the corresponding entry from P and start reverse bumping
            x = P[i].pop(j)
            if len(P[i]) == 0:
                P.pop(i)

            # move upwards: in row r find rightmost entry < x, replace it and continue,
            # if none found in a row, stop and output x
            for r in range(i - 1, -1, -1):
                replaced = False
                for k in range(len(P[r]) - 1, -1, -1):
                    if P[r][k] < x:
                        y = P[r][k]
                        P[r][k] = x
                        x = y
                        replaced = True
                        break
                if not replaced:
                    break

            word_rev.append(int(x))

        # reverse collected letters to obtain original insertion order
        return list(reversed(word_rev))

    @classmethod
    def rsk_insert(cls, *letters):
        """
        Perform ordinary RSK (row insertion) on the given sequence of letters,
        starting from this Plactic as the initial P-tableau. Returns a pair
        (P_tableau, Q_tableau) where both are Plactic instances and Q is the
        standard recording tableau with entries 1..m (in insertion order).

        Usage:
          P, Q = Plactic().rsk_insert(3,1,2,1)
          or
          P, Q = Plactic().rsk_insert([3,1,2,1])
        """
        # normalize letters accept either rsk_insert(a,b,...) or rsk_insert(iterable)
        if len(letters) == 1 and isinstance(letters[0], (list, tuple)):
            seq = list(letters[0])
        else:
            seq = list(letters)

        # mutable rows for P and Q
        P_rows = []
        Q_rows = []

        def _insert_letter(rows, x):
            """Insert x into mutable rows (list of list) and return final (i,j)."""
            i = 0
            while True:
                if i >= len(rows):
                    rows.append([])
                row = rows[i]
                if len(row) == 0 or x >= max(row):
                    row.append(x)
                    return i, len(row) - 1
                # bump smallest entry > x
                for idx, val in enumerate(row):
                    if val > x:
                        bumped = val
                        row[idx] = x
                        x = bumped
                        break
                i += 1

        # perform insertions, updating Q at appended cell each time
        for k, letter in enumerate(seq, start=1):
            i, j = _insert_letter(P_rows, int(letter))
            # ensure Q has row i
            while i >= len(Q_rows):
                Q_rows.append([])
            # append placeholder zeros for earlier existing length mismatches
            # (Q_rows[i] should be same length as P_rows[i]-1 before append)
            # append the insertion index to the new/extended cell
            if j == len(Q_rows[i]):
                Q_rows[i].append(k)
            else:
                # In typical RSK we always append a new box at the final row.
                # However guard against unexpected states by inserting at position j.
                Q_rows[i].insert(j, k)
                # ensure consistency of row lengths in Q with P:
                # if this caused earlier rows to mismatch, leave as-is (shouldn't happen)
        # convert to Plactic objects
        P = cls(tuple(tuple(r) for r in P_rows))
        Q = cls(tuple(tuple(r) for r in Q_rows))
        return P, Q

    def reverse_rectify_to_outer(self, outer_shape):
        r"""
        Deterministic reverse-rectification to a given outer shape `outer_shape`.

        Given a (straight) tableau `self` of shape lambda, produce a skew tableau
        (represented as a Plactic whose rows may contain 0's for inner cells)
        of outer shape `outer_shape` whose rectification is `self`.

        outer_shape: iterable of nonnegative ints giving the desired outer row
                     lengths (mu_0 >= mu_1 >= ...).

        Algorithm (deterministic):
        - Let mu be the set of cells (r,c) with 0 <= r < len(mu) and 0 <= c < mu[r].
        - While the current set of occupied cells (from the working tableau) is
          a strict subset of mu:
            * choose an outer corner cell (r,c) in mu\current_cells (no cell of mu
              to its right or below). Choose the maximal such (r,c) (deterministic).
            * create a hole at (r,c) (extend rows/cols as needed, set that cell to 0),
              then perform an upward jeu-de-taquin slide from (r,c) using
              up_jdt_slide to move the hole inward.
            * adopt the resulting tableau and continue.
        - Return the resulting Plactic (with zeros marking inner/removed cells).

        Notes:
        - Raises ValueError if outer_shape does not dominate the current shape
          (i.e. mu must contain the current occupied cells).
        - Raises RuntimeError if no suitable outer corner can be found or if an
          up_jdt_slide fails (this indicates the requested outer shape is not
          attainable by reverse-rectification).
        """
        # normalize outer shape
        mu = tuple(int(x) for x in outer_shape)
        if any(mu[i] < (mu[i + 1] if i + 1 < len(mu) else 0) for i in range(len(mu))):
            # not required but warn if not nonincreasing
            pass

        # Convert grid to mutable list of lists (excluding border)
        rows = []
        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                val = self._grid[i, j]
                if val is not None:
                    row.append(val if val != 0 else 0)
            if row or i < len(mu):  # Keep rows that might be needed
                rows.append(row)

        # helper to compute occupied cells set
        def occupied_cells(rs):
            cells = set()
            for rr, row in enumerate(rs):
                for cc in range(len(row)):
                    # treat only real boxes (non-zero and zero both count as occupied for shape)
                    cells.add((rr, cc))
            return cells

        # build mu cell set
        mu_cells = set()
        for r in range(len(mu)):
            for c in range(mu[r]):
                mu_cells.add((r, c))

        cur_cells = occupied_cells(rows)

        # sanity: current cells must be subset of mu_cells
        if not cur_cells.issubset(mu_cells):
            raise ValueError("outer_shape must contain the current tableau shape (componentwise domination)")

        # loop until shapes match
        while cur_cells != mu_cells:
            # difference cells that need to be added
            diff = mu_cells - cur_cells

            # find outer-corner candidates inside diff (no mu cell to right or below)
            candidates = []
            for r, c in diff:
                right = (r, c + 1) in mu_cells
                down = (r + 1, c) in mu_cells
                if (not right) and (not down):
                    candidates.append((r, c))

            if not candidates:
                raise RuntimeError("No outer-corner candidate found while reverse-rectifying to outer shape")

            # deterministic choice: maximal (r,c)
            r, c = max(candidates)

            # ensure rows has row r
            while len(rows) <= r:
                rows.append([])

            # extend the chosen row to have column c (fill with zeros)
            while len(rows[r]) <= c:
                rows[r].append(0)

            # build temporary tableau with the new hole and perform up_jdt_slide
            tmp = Plactic(tuple(tuple(rr) for rr in rows))

            # verify starting position is empty (should be 0)
            val = tmp[r, c]
            if val is None:
                # if indexing returned None, ensure underlying grid has that slot
                # Set to 0 in the working rows
                rows[r][c] = 0
                val = 0
            if val != 0 and val is not None:
                # force the hole (overwrite); we must set it to 0 before sliding
                rows[r][c] = 0
                tmp = Plactic(tuple(tuple(rr) for rr in rows))

            # perform upward slide from the hole
            try:
                new_tmp = tmp.up_jdt_slide(r, c)
            except Exception as e:
                raise RuntimeError(f"up_jdt_slide failed from corner {(r, c)}: {e}")

            # adopt new rows and recompute occupied cells
            # Convert grid back to list of lists (excluding border)
            rows = []
            for i in range(new_tmp.rows):
                row = []
                for j in range(new_tmp.cols):
                    val = new_tmp._grid[i, j]
                    if val is not None:
                        row.append(val if val != 0 else 0)
                if row or i < len(mu):
                    rows.append(row)
            cur_cells = occupied_cells(rows)

        # return Plactic with zeros possibly marking inner cells
        return Plactic(tuple(tuple(r) for r in rows))

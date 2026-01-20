"""
Bumpless Pipe Dreams (BPD) module
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import IntEnum
from functools import cache, cached_property
from typing import Tuple

import numpy as np
from sympy import pretty
from sympy.printing.defaults import DefaultPrinting

from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.schub_lib.schubert_monomial_graph import SchubertMonomialGraph
from schubmult.symbolic import Expr


class TileType(IntEnum):
    """
    Enumeration of the 6 possible tile types in a pipe dream.

    Each tile represents how two pipes (horizontal and vertical) interact in a square.
    """

    TBD = 0  # Placeholder for uninitialized tile
    BLANK = 1  # Both pipes go straight (no crossing, no elbow)
    CROSS = 2  # Pipes cross each other
    HORIZ = 3
    ELBOW_NW = 4  # Elbow: bottom-right to top-left (╯)
    ELBOW_SE = 5  # Elbow: top-left to bottom-right (╮)
    VERT = 6
    BUMP = 7  # Bump/osculating tile (pipes touch at corner)

    def __str__(self) -> str:
        symbols = {TileType.BLANK: "▢", TileType.CROSS: "┼", TileType.ELBOW_NW: "╯", TileType.ELBOW_SE: "╭", TileType.HORIZ: "─", TileType.VERT: "│", TileType.BUMP: "╬", TileType.TBD: "?"}
        return symbols.get(self, "?")

    @cached_property
    def is_crossing(self) -> bool:
        """True if this tile is a crossing"""
        return self == TileType.CROSS

    @cached_property
    def is_elbow(self) -> bool:
        """True if this tile is any type of elbow"""
        return self in (TileType.ELBOW_NW, TileType.ELBOW_SE)

    @cached_property
    def is_empty(self) -> bool:
        """True if this tile is empty (pipes go straight)"""
        return self == TileType.BLANK

    @cached_property
    def feeds_right(self) -> bool:
        """True if the horizontal pipe continues to the right"""
        return self in (TileType.HORIZ, TileType.ELBOW_SE, TileType.CROSS, TileType.BUMP)

    @cached_property
    def feeds_up(self) -> bool:
        """True if the vertical pipe continues upwards"""
        return self in (TileType.VERT, TileType.ELBOW_NW, TileType.CROSS, TileType.BUMP)

    @cached_property
    def entrance_from_bottom(self) -> bool:
        """True if a pipe can enter from the bottom"""
        return self in (TileType.VERT, TileType.ELBOW_SE, TileType.CROSS, TileType.BUMP)

    @cached_property
    def entrance_from_left(self) -> bool:
        """True if a pipe can enter from the left"""
        return self in (TileType.HORIZ, TileType.ELBOW_NW, TileType.CROSS, TileType.BUMP)


def _invalidate_grid(grid: np.ndarray) -> None:
    """Helper function to invalidate a grid by setting TBD tiles"""
    # Compute mask once: keep only BLANK and CROSS tiles
    mask = (grid == TileType.BLANK) | (grid == TileType.CROSS) | (grid == TileType.BUMP)
    grid[~mask] = TileType.TBD


# def _is_asm(arr):
#     # Ensure it's a square numpy array
#     n, m = arr.shape
#     if n != m:
#         return False

#     # 1. Check Row and Column Sums (Must be exactly 1)
#     if not np.all(arr.sum(axis=0) == 1) or not np.all(arr.sum(axis=1) == 1):
#         return False

#     # 2. Check Allowed Values (Only 0, 1, -1)
#     # Using absolute value is a fast way to check this
#     # if not np.all(np.abs(arr) <= 1):
#     #     return False

#     # 3. Check Alternating Signs and Boundary +1s
#     # In a valid ASM, the cumulative sum along any row/column
#     # must stay between 0 and 1 at all times.
#     # Because row sums are 1, the first/last non-zero must be 1
#     # if the prefix sums are always 0 or 1.

#     # Check rows
#     row_cumsum = np.cumsum(arr, axis=1)
#     if np.any((row_cumsum < 0) | (row_cumsum > 1)):
#         return False

#     # Check columns
#     col_cumsum = np.cumsum(arr, axis=0)
#     if np.any((col_cumsum < 0) | (col_cumsum > 1)):
#         return False

#     # 4. Final verification: The non-zero entries must actually change.
#     # Example of a failure mode for cumsum alone: [1, 0, 0] is fine,
#     # but [1, 1, -1] has cumsums [1, 2, 1], which is caught above.
#     # The only thing left is to ensure we don't have two 1s or two -1s in a row.
#     # We can do this by checking if the absolute difference of non-zero entries is 2.


#     # Vectorized adjacent nonzero check for rows and columns
#     for axis in [0, 1]:
#         lines = arr if axis == 1 else arr.T
#         # Create a mask of nonzero entries
#         mask = lines != 0
#         # Pad with False at the end to align for diff
#         mask_shifted = np.roll(mask, -1, axis=1)
#         mask_shifted[:, -1] = False
#         # Find indices of nonzero entries
#         # For each line, get the nonzero values and check adjacent diffs
#         # We'll use a trick: for each line, get the indices of nonzeros, then check if any adjacent values are equal
#         # This is still a little tricky to fully vectorize, but we can batch it with list comprehension
#         if np.any([
#             np.any(np.diff(line[line != 0]) == 0)
#             for line in lines
#         ]):
#             return False

#     return True


class BPD(SchubertMonomialGraph, DefaultPrinting):
    """
    Bumpless Pipe Dream representation.

    A bumpless pipe dream is an n×n grid where:
    - TileType.CROSS (1) represents a crossing
    - TileType.BLANK (0) represents an empty box (pipes go straight)
    - For general pipe dreams, can use TileType.ELBOW_* (2-5) for elbows

    Each BPD corresponds to a permutation and has an associated weight.
    """

    def __init__(self, grid, column_perm: Permutation | None = None, *, _is_copy=False) -> None:
        """
        Initialize a BPD from a grid.

        Args:
            grid: n×n array-like of TileType values, integers 0-5, or list of lists
        s"""
        if _is_copy:
            return
        self._grid = np.array(grid, dtype=TileType)
        self._column_perm = column_perm if column_perm else Permutation([])
        self._perm = None
        self._valid = None
        self._word = None
        self._unzero_cache = None
        self.build()

    @property
    def rows(self) -> int:
        return self._grid.shape[0]

    @property
    def cols(self) -> int:
        return self._grid.shape[1]

    @classmethod
    @cache
    def all_bpds(cls, w: Permutation, length: int | None = None) -> set[BPD]:
        if length is None:
            length = len(w)
        pipes = set()
        new_pipes = [BPD.rothe_bpd(w, length)]

        while len(new_pipes) != 0:
            pipe = new_pipes.pop(0)
            pipes.add(pipe)

            for move in pipe.droop_moves():
                new_pipe = pipe.do_droop_move(move)
                if new_pipe not in pipes:
                    new_pipes.append(new_pipe)
                    pipes.add(new_pipe)
        return pipes

    # Static lookup table for TBD tile resolution
    # Index by (left_tile_value, up_tile_value) where None is represented as -1
    # Pre-computed based on feeds_right and entrance_from_bottom properties
    _TBD_LOOKUP = None

    @classmethod
    def _init_tbd_lookup(cls):
        """Initialize the static TBD tile lookup table."""
        if cls._TBD_LOOKUP is not None:
            return

        # Create lookup table: shape (9, 9) for tile types 0-7 (including BUMP=7) plus None (-1 -> index 8)
        lookup = np.full((9, 9), TileType.ELBOW_SE, dtype=TileType)

        # Map None to index 8
        none_idx = 8

        # For each combination of (left_tile, up_tile), compute result
        # Range includes -1 (None) through 7 (BUMP)
        for left_val in range(-1, 8):  # -1 for None, 0-7 for tile types (TBD, BLANK, CROSS, HORIZ, ELBOW_NW, ELBOW_SE, VERT, BUMP)
            left_idx = none_idx if left_val == -1 else left_val
            left_feeds_right = False if left_val == -1 else TileType(left_val).feeds_right

            for up_val in range(-1, 8):
                up_idx = none_idx if up_val == -1 else up_val
                up_entrance = False if up_val == -1 else TileType(up_val).entrance_from_bottom

                # Apply the logic from _get_tbd_tile
                if up_val == -1:  # up_tile is None
                    if left_val == -1 or not left_feeds_right:
                        lookup[left_idx, up_idx] = TileType.ELBOW_SE
                    else:
                        lookup[left_idx, up_idx] = TileType.HORIZ
                elif left_val == -1:  # left_tile is None, up_tile is not
                    if up_entrance:
                        lookup[left_idx, up_idx] = TileType.VERT
                    else:
                        lookup[left_idx, up_idx] = TileType.ELBOW_SE
                else:  # Both not None
                    if left_feeds_right:
                        if up_entrance:
                            lookup[left_idx, up_idx] = TileType.ELBOW_NW
                        else:
                            lookup[left_idx, up_idx] = TileType.HORIZ
                    else:
                        if up_entrance:
                            lookup[left_idx, up_idx] = TileType.VERT
                        else:
                            lookup[left_idx, up_idx] = TileType.ELBOW_SE

        cls._TBD_LOOKUP = lookup

    @classmethod
    def _get_tbd_tile(cls, left_tile: TileType | None, up_tile: TileType | None) -> TileType:
        """Fast lookup-based TBD tile resolution."""
        if cls._TBD_LOOKUP is None:
            cls._init_tbd_lookup()

        left_idx = 8 if left_tile is None else int(left_tile)
        up_idx = 8 if up_tile is None else int(up_tile)
        return cls._TBD_LOOKUP[left_idx, up_idx]

    DEBUG = False

    def build(self) -> None:
        """Build internal structures by resolving TBD tiles using lookup table."""
        # Cache shape to avoid property access overhead
        rows = self._grid.shape[0]
        cols = self._grid.shape[1]

        if rows == 0:
            return

        # Initialize lookup table if needed
        if BPD._TBD_LOOKUP is None:
            BPD._init_tbd_lookup()

        # Find all TBD positions
        tbd_mask = self._grid == TileType.TBD

        # Corner case [0,0]
        if tbd_mask[0, 0]:
            self._grid[0, 0] = BPD._TBD_LOOKUP[8, 8]  # (None, None)

        # First column [1:, 0] - must process sequentially as each row depends on previous
        if rows > 1:
            tbd_rows = np.where(tbd_mask[1:, 0])[0] + 1
            for row in tbd_rows:
                up_tile = int(self._grid[row - 1, 0])
                self._grid[row, 0] = BPD._TBD_LOOKUP[8, up_tile]

        # Process column by column (dependencies require sequential processing)
        for col in range(1, cols):
            # Top row [0, col] - single lookup with (left_tile, None)
            if tbd_mask[0, col]:
                left_tile = int(self._grid[0, col - 1])
                self._grid[0, col] = BPD._TBD_LOOKUP[left_tile, 8]

            # Interior cells [1:, col] - must process sequentially as each row depends on previous
            if rows > 1:
                tbd_rows = np.where(tbd_mask[1:, col])[0] + 1
                for row in tbd_rows:
                    left_tile = int(self._grid[row, col - 1])
                    up_tile = int(self._grid[row - 1, col])
                    self._grid[row, col] = BPD._TBD_LOOKUP[left_tile, up_tile]

    def __len__(self) -> int:
        """Return the size n of the n×n grid"""
        return self.rows

    def __getitem__(self, key) -> TileType | np.ndarray:
        """Access grid elements, casting to TileType"""
        return self._grid[key]

    def _sympyrepr(self, printer=None) -> str:
        """SymPy repr representation"""
        grid_list = self._grid.tolist()
        if printer is None:
            return f"BPD({grid_list})"
        return f"BPD({printer._print(grid_list)})"

    def _sympystr(self, printer=None) -> str:
        """SymPy str representation of the BPD using tile symbols"""
        return printer._print("BPD(\n" + pretty(self) + ")")

    def _pretty(self, printer=None):
        """Pretty printing with row and column labels"""
        from sympy.printing.pretty.stringpict import prettyForm

        if printer is None:
            # Fallback to simple string representation
            return prettyForm(self._sympystr(None))

        # Calculate column widths based on labels (minimum 1 for tile)
        col_widths = []
        for j in range(self.cols):
            label = str(printer._print(self._column_perm[j]))
            col_widths.append(max(1, len(label)))

        # Use maximum column width for all columns
        max_col_width = max(col_widths)

        # Build rows with proper tile extensions
        rows = []
        perm = self.perm
        perm_values = [perm[i] for i in range(self.rows)]

        for i in range(self.rows):
            row_parts = []
            for j in range(self.cols):
                tile = self[i, j]
                tile_str = str(tile)

                # Symmetric padding to center tile in column using max width
                total_pad = max(max_col_width - 1, 2)
                left_pad = total_pad // 2
                right_pad = total_pad - left_pad

                horiz_str = str(TileType.HORIZ)
                space_str = " "

                # Determine left padding character
                if tile.entrance_from_left:
                    left_str = horiz_str * left_pad
                else:
                    left_str = space_str * left_pad

                # Determine right padding character
                if tile.feeds_right:
                    right_str = horiz_str * right_pad
                else:
                    right_str = space_str * right_pad

                tile_str = left_str + tile_str + right_str

                row_parts.append(tile_str)

            row_str = "".join(row_parts)
            # Add permutation value as row label on the right
            perm_label = str(printer._print(perm_values[i]))
            row_str += f" {perm_label}"
            rows.append(row_str)

        # Add column labels at the bottom with proper spacing and alignment
        col_labels = []
        # Each column is 1 (tile) + total_pad wide
        total_pad = max(max_col_width - 1, 2)
        column_width = 1 + total_pad

        for j in range(self.cols):
            label = str(printer._print(self._column_perm[j]))
            col_labels.append(label.center(column_width))

        # Add leading space to match prettyForm's handling
        rows.append("".join(col_labels) + (" " * max_col_width) + " ")

        return prettyForm("\n".join(rows))

    def shiftup(self, shift: int = 1) -> BPD:
        """Shift the BPD up by a given amount."""
        # Create new grid with shifted dimensions
        new_rows = self.rows + shift
        new_cols = self.cols + shift
        new_grid = np.full((new_rows, new_cols), TileType.ELBOW_SE, dtype=TileType)

        # Shift the grid contents
        for i in range(self.rows):
            for j in range(self.cols):
                new_grid[i + shift, j + shift] = self[i, j]

        # Fill the top-left portion with identity pattern
        for i in range(shift):
            for j in range(shift):
                if i == j:
                    new_grid[i, j] = TileType.ELBOW_SE
                elif i < j:
                    new_grid[i, j] = TileType.HORIZ
                else:
                    new_grid[i, j] = TileType.VERT

        # Connect the identity to shifted content
        for i in range(shift):
            for j in range(shift, new_cols):
                new_grid[i, j] = TileType.HORIZ
        for i in range(shift, new_rows):
            for j in range(shift):
                new_grid[i, j] = TileType.VERT
        return BPD(new_grid)

    @property
    def perm(self) -> Permutation:
        """
        Compute the permutation associated with this BPD.

        The permutation is determined by following each vertical pipe from bottom to top.
        Pipes enter from the bottom (vertical) and left (horizontal).

        Returns:
            Permutation object
        """
        if self._perm is not None:
            return self._perm
        if self.rows == 0:
            self._perm = Permutation([])
            return self._perm
        # self._perm = Permutation.ref_product(*self.word)

        nrows, ncols = self._grid.shape
        bottom_row = self._grid[nrows - 1, :]

        # Check for TBD in bottom row
        if np.any(bottom_row == TileType.TBD):
            raise ValueError("Cannot compute permutation with unresolved TBD tiles")

        # Vectorized: find columns with entrance_from_bottom
        # entrance_from_bottom is True for VERT, CROSS, ELBOW_NW
        entrance_mask = (bottom_row == TileType.VERT) | (bottom_row == TileType.CROSS) | (bottom_row == TileType.ELBOW_SE)
        good_cols = (np.where(entrance_mask)[0] + 1).tolist()

        good_cols = Permutation.from_partial(good_cols)

        small_perm = Permutation([])
        # Vectorized: Map tiles to their diff values
        diff = np.ones((nrows, ncols), dtype=int)
        diff[self._grid == TileType.BLANK] = 0
        diff[self._grid == TileType.CROSS] = 2
        # Create r array with shape (nrows+1, ncols+1)
        r = np.zeros((nrows + 1, ncols + 1), dtype=int)

        for i in range(1, nrows + 1):
            for j in range(1, ncols + 1):
                r[i, j] = r[i - 1, j - 1] + diff[i - 1, j - 1]

        # Pre-compute all cross positions and their pipes_northeast values
        cross_positions = np.argwhere(self._grid == TileType.CROSS)
        if len(cross_positions) > 0:
            # Sort by column first, then by row descending (for correct swap order)
            sort_indices = np.lexsort((-cross_positions[:, 0], cross_positions[:, 1]))
            cross_positions = cross_positions[sort_indices]
            # Get pipes_northeast for each cross position
            pipes_northeast_values = r[cross_positions[:, 0] + 1, cross_positions[:, 1] + 1]
            # Apply swaps sequentially
            for pipes_northeast in pipes_northeast_values:
                small_perm = small_perm.swap(pipes_northeast - 2, pipes_northeast - 1)

        build_perm = good_cols * small_perm

        self._perm = Permutation.from_partial(build_perm)
        return self._perm

    @property
    def permutation(self) -> Permutation:
        """Alias for perm property"""
        return self.perm

    @property
    def inv(self) -> int:
        """
        Return the inversion count of the associated permutation.

        This is a convenience property that delegates to perm.inv.

        Returns:
            Number of inversions in the permutation
        """
        return self.perm.inv

    @property
    def length_vector(self) -> tuple[int, ...]:
        """
        Compute the length vector of the permutation represented by this BPD.

        The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
        of crossings in row i.

        Returns:
            Tuple of integers representing the length vector
        """
        # Vectorized: compute all row sums at once
        blank_counts = np.sum(self._grid == TileType.BLANK, axis=1)
        return tuple(int(x) for x in blank_counts)

    @classmethod
    def from_asm(cls, asm) -> BPD:
        """
        Create a BPD from an ASM (Alternating Sign Matrix).

        Args:
            asm: n×n array-like of integers (-1, 0, 1)
        Returns:
            BPD object
        """
        corner_sum = np.cumsum(np.cumsum(asm, axis=0), axis=1)
        n = asm.shape[0]
        grid = np.full((n, n), fill_value=TileType.TBD, dtype=TileType)

        def asm_tile(i, j):
            if corner_sum[i, j] - (corner_sum[i - 1, j - 1] if i > 0 and j > 0 else 0) == 2:
                return TileType.CROSS
            if corner_sum[i, j] - (corner_sum[i - 1, j - 1] if i > 0 and j > 0 else 0) == 0:
                return TileType.BLANK
            if (corner_sum[i - 1, j] if i > 0 else 0) - (corner_sum[i, j - 1] if j > 0 else 0) == -1:
                return TileType.HORIZ
            if (corner_sum[i - 1, j] if i > 0 else 0) - (corner_sum[i, j - 1] if j > 0 else 0) == 1:
                return TileType.VERT

            raise ValueError(f"Invalid ASM corner sum at ({i}, {j})")

        asm = np.array(asm, dtype=np.int8)
        grid[asm == 1] = TileType.ELBOW_SE
        grid[asm == -1] = TileType.ELBOW_NW
        rows, cols = np.where(asm == 0)
        grid[rows, cols] = np.vectorize(asm_tile, otypes=[TileType])(rows, cols)
        return cls(grid)

    def to_asm(self):
        to_trans = self.normalize()
        n = max(to_trans.rows, to_trans.cols)
        asm = np.full(shape=(n, n), fill_value=0, dtype=np.int8)
        se_elbows = to_trans.all_se_elbows()
        if len(se_elbows) > 0:
            rows, cols = zip(*se_elbows)
            asm[rows, cols] = 1
        nw_elbows = to_trans.all_nw_elbows()
        if len(nw_elbows) > 0:
            rows, cols = zip(*nw_elbows)
            asm[rows, cols] = -1
        # print(asm)
        return asm

    @classmethod
    @cache
    def rothe_bpd(cls, perm: Permutation, num_rows: int | None = None) -> BPD:
        if num_rows is None:
            num_rows = len(perm)
        n = max(num_rows, len(perm))
        grid = np.full((num_rows, n), fill_value=TileType.TBD, dtype=TileType)

        # Set blanks from diagram
        diagram = perm.diagram
        if diagram:
            diagram_arr = np.array(list(diagram), dtype=np.int32)
            grid[diagram_arr[:, 0] - 1, diagram_arr[:, 1] - 1] = TileType.BLANK

        # Vectorize cross detection using graph coordinates
        graph = perm.graph
        if graph and len(graph) > 0:
            graph_arr = np.array(list(graph), dtype=np.int32)

            # Get mask of non-blank positions
            non_blank_mask = grid != TileType.BLANK
            non_blank_coords = np.argwhere(non_blank_mask)

            if len(non_blank_coords) > 0:
                # For each non-blank position, check graph conditions
                for i, j in non_blank_coords:
                    # Check if any graph point has row == i+1 and col < j+1
                    has_row_below = np.any((graph_arr[:, 0] == i + 1) & (graph_arr[:, 1] - 1 < j))
                    # Check if any graph point has row < i+1 and col == j+1
                    has_col_left = np.any((graph_arr[:, 0] - 1 < i) & (graph_arr[:, 1] == j + 1))
                    if has_row_below and has_col_left:
                        grid[i, j] = TileType.CROSS

        return BPD(grid)

    @property
    def weight(self) -> Tuple[int, ...]:
        """
        Compute the weight of this BPD.

        The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
        of empty squares (0s) in column i.

        Returns:
            Tuple of integers representing the weight
        """
        return tuple(int(np.sum(self[:, j] == TileType.BLANK)) for j in range(self.rows))

    @property
    def word(self) -> Tuple[int, ...]:
        """
        Compute a reduced word for the permutation represented by this BPD.

        For each crossing at position (i,j), the word value is the number of pipes
        weakly northeast of the crossing minus 1. Weakly northeast means all positions
        (r,c) where r <= i and c >= j.

        Returns:
            Tuple of integers representing the reduced word (1-indexed positions)
        """
        if self._word is not None:
            return self._word

        word = []

        nrows, ncols = self._grid.shape
        # Map tiles to their diff values
        diff = np.ones_like(self._grid, dtype=int)
        diff[self._grid == TileType.BLANK] = 0
        diff[self._grid == TileType.CROSS] = 2
        # Create r array with shape (nrows+1, ncols+1)
        r = np.zeros((nrows + 1, ncols + 1), dtype=int)

        for i in range(1, nrows + 1):
            for j in range(1, ncols + 1):
                r[i, j] = r[i - 1, j - 1] + diff[i - 1, j - 1]
        # return r[target_row + 1, target_col + 1]
        for col in range(self.cols):
            for row in range(self.rows - 1, -1, -1):
                if self[row, col] == TileType.CROSS:
                    pipes_northeast = r[row + 1, col + 1]  # self.cols] - r[row + 1, col]
                    word.append(pipes_northeast - 1)
        self._word = tuple(word)
        return self._word

    def set_width(self, width):
        """Set the width of the BPD by adding empty columns on the right if needed."""
        if width < self.cols:
            if len(self.perm) > width:
                raise ValueError("New width must be at least the length of the permutation")
            new_grid = np.full((self.rows, width), fill_value=TileType.TBD, dtype=TileType)
            new_grid[:, :width] = self._grid[:, :width]
            bop = BPD(new_grid)
            bop.rebuild()
            return bop
        if width == self.cols:
            return self
        new_grid = np.full((self.rows, width), fill_value=TileType.TBD, dtype=TileType)
        new_grid[:, : self.cols] = self._grid
        return BPD(new_grid)

    @property
    def is_valid(self) -> bool:
        """
        Check if this is a valid bpd.

        Returns:
            True if valid, False otherwise
        """
        if self._valid is not None:
            return self._valid
        if self.rows == 0:
            self._valid = True
            return self._valid
        # if self.perm.inv != sum(self.length_vector):
        #     self._valid = False
        #     return self._valid
        if any(self[0, j].feeds_up for j in range(self.cols)):
            self._valid = False
            return self._valid
        if any(self[i, 0].entrance_from_left for i in range(self.rows)):
            self._valid = False
            return self._valid
        if any(not self[i, self.cols - 1].feeds_right for i in range(self.rows)):
            self._valid = False
            return self._valid
        for i in range(1, self.rows):
            for j in range(1, self.cols - 1):
                if self[i, j].entrance_from_left and not self[i, j - 1].feeds_right:
                    self._valid = False
                    return self._valid
                if self[i, j].feeds_up and not self[i - 1, j].entrance_from_bottom:
                    self._valid = False
                    return self._valid
                if self[i, j].feeds_right and not self[i, j + 1].entrance_from_left:
                    self._valid = False
                    return self._valid
                if self[i, j].entrance_from_bottom and i < self.rows - 1 and not self[i + 1, j].feeds_up:
                    self._valid = False
                    return self._valid
        if len(self.perm.trimcode) > self.rows:
            self._valid = False
            return self._valid
        self._valid = True
        return self._valid

    def __eq__(self, other: object) -> bool:
        """Check equality of two BPDs"""
        if not isinstance(other, BPD):
            return False
        return np.array_equal(self._grid, other._grid)

    def __hash__(self) -> int:
        """Hash for use in sets and dicts"""
        return hash(self._grid.tobytes())

    def copy(self) -> BPD:
        """Create a copy of this BPD"""
        new_bpd = BPD(None, _is_copy=True)
        new_bpd._grid = self._grid.copy()
        new_bpd._column_perm = self._column_perm
        new_bpd._perm = self._perm
        new_bpd._valid = self._valid
        new_bpd._word = self._word
        return new_bpd

    @property
    def num_crossings(self) -> int:
        """Total number of crossings in the BPD"""
        return int(np.sum(self._grid == TileType.CROSS))

    def inversion_at_cross(self, i: int, j: int) -> int:
        """
        Compute the inversion associated with the crossing at position (i, j).

        The inversion is determined by tracing the pipes through the BPD.

        Args:
            i: Row index of the crossing
            j: Column index of the crossing
        Returns:
            The inversion value as an integer
        """
        if self[i, j] != TileType.CROSS:
            return None

        # up and right
        up_r = i - 1
        up_c = j

        while up_r >= 0 and up_c < self.cols:
            tile = self[up_r, up_c]
            if tile.feeds_up:
                up_r -= 1
            else:
                up_c += 1

        right_r = i
        right_c = j + 1

        while right_r >= 0 and right_c < self.cols:
            tile = self[right_r, right_c]
            if tile.feeds_right:
                right_c += 1
            else:
                right_r -= 1
        assert up_r < right_r
        return up_r, right_r

    def trace_pipe(self, i: int, j: int, direction: str | None = None) -> int | None:
        if self[i, j] == TileType.ELBOW_NW:
            if direction == "left":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self[i, j] == TileType.ELBOW_SE:
            if direction == "down":
                return None
            if i == self.rows - 1:
                return self._column_perm[j]
            return self.trace_pipe(i + 1, j, direction="down")
        if self[i, j] == TileType.HORIZ:
            if direction == "down":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self[i, j] == TileType.VERT:
            if direction == "left":
                return None
            if i == self.rows - 1:
                return self._column_perm[j]
            return self.trace_pipe(i + 1, j, direction="down")
        if self[i, j] == TileType.CROSS:
            if direction == "down":
                if i == self.rows - 1:
                    return self._column_perm[j]
                return self.trace_pipe(i + 1, j, direction="down")
            if direction == "left":
                return self.trace_pipe(i, j - 1, direction="left")
            raise ValueError("Must specify direction when tracing through a crossing")
        if self[i, j] == TileType.BUMP:
            if direction == "down":
                return self.trace_pipe(i, j - 1, direction="left")
            if direction == "left":
                if i == self.rows - 1:
                    return self._column_perm[j]
                return self.trace_pipe(i + 1, j, direction="down")
            raise ValueError("Must specify direction when tracing through a crossing")
        raise ValueError(f"Invalid tile for tracing pipe at ({i}, {j}): {self[i, j]}")

    def all_se_elbows(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.ELBOW_SE)

    def all_nw_elbows(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.ELBOW_NW)

    def all_blanks(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.BLANK)

    def all_crossings(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.CROSS)

    def all_tiles_of_type(self, tile_type: TileType) -> set[tuple[int, int]]:
        if isinstance(tile_type, (list, tuple)):
            result = set()
            for t in tile_type:
                result.update(self.all_tiles_of_type(t))
            return result
        return set({(int(i), int(j)) for i, j in zip(*np.where(self._grid == tile_type))})

    def droop_moves(self) -> set[tuple[tuple[int, int], tuple[int, int]]]:
        import itertools

        droop_moves = set()
        for (ri, rj), (bi, bj) in itertools.product(self.all_se_elbows(), self.all_blanks()):
            if bi > ri and bj > rj:
                LEGAL = True
                for i, j in itertools.product(range(ri, bi + 1), range(rj, bj + 1)):
                    if (i, j) != (ri, rj) and (self[i, j] == TileType.ELBOW_SE or self[i, j] == TileType.ELBOW_NW or self[i, j] == TileType.BUMP):  # if not NW-corner, check if elbow
                        LEGAL = False
                        break
                if LEGAL:
                    droop_moves.add(((ri, rj), (bi, bj)))
        return droop_moves

    def min_droop_moves(self) -> set[tuple[tuple[int, int], tuple[int, int]]]:
        droop_moves = set()
        ARBITRARY_LARGE = 200
        for a, b in self.all_tiles_of_type((TileType.ELBOW_SE, TileType.BUMP)):
            x = np.min(np.argwhere(self._grid[a + 1 :, b] != TileType.CROSS), initial=ARBITRARY_LARGE) + a + 1
            if x > ARBITRARY_LARGE:
                continue
            y = np.min(np.argwhere(self._grid[a, b + 1 :] != TileType.CROSS), initial=ARBITRARY_LARGE) + b + 1
            if y > ARBITRARY_LARGE:
                continue
            droop_moves.add(((a, b), (x, y)))
        return droop_moves

    def do_min_droop_move(self, move: tuple[tuple[int, int], tuple[int, int]]) -> BPD:
        D = self.copy()
        (ri, rj) = move[0]
        (bi, bj) = move[1]

        # if bi >= D.rows:
        #     D = D.resize(bi + 2, check_valid=False)
        # if bj >= D.cols:
        #     D = D.set_width(bj + 2)
        orig_self = self.resize(D.rows, check_valid=False).set_width(D.cols)
        D._grid[ri, rj] = TileType.BLANK if orig_self[ri, rj] == TileType.ELBOW_SE else TileType.ELBOW_NW  # NW-corner
        D._grid[bi, bj] = TileType.ELBOW_NW if orig_self[bi, bj] == TileType.BLANK else TileType.BUMP  # SE-corner
        D._grid[bi, rj] = TileType.ELBOW_SE  # SW-corner
        D._grid[ri, bj] = TileType.ELBOW_SE  # NE-corner

        # top and bottom
        for j in range(rj + 1, bj):
            if orig_self[ri, j] == TileType.HORIZ:
                D._grid[ri, j] = TileType.BLANK
            else:  # self[ri, j] == TileType.CROSS
                D._grid[ri, j] = TileType.VERT
            if orig_self[bi, j] == TileType.BLANK:
                D._grid[bi, j] = TileType.HORIZ
            else:  # self[bi, j] == TileType.VERT
                D._grid[bi, j] = TileType.CROSS
        # left and right
        for i in range(ri + 1, bi):
            if orig_self[i, rj] == TileType.VERT:
                D._grid[i, rj] = TileType.BLANK
            else:  # self[i, rj] == TileType.CROSS
                D._grid[i, rj] = TileType.HORIZ
            if orig_self[i, bj] == TileType.BLANK:
                D._grid[i, bj] = TileType.VERT
            else:
                D._grid[i, bj] = TileType.CROSS
        return D

    def do_droop_move(self, move: tuple[tuple[int, int], tuple[int, int]]) -> BPD:
        D = self.copy()
        (ri, rj) = move[0]
        (bi, bj) = move[1]

        D._grid[ri, rj] = TileType.BLANK  # NW-corner
        D._grid[bi, bj] = TileType.ELBOW_NW  # SE-corner
        D._grid[bi, rj] = TileType.ELBOW_SE  # SW-corner
        D._grid[ri, bj] = TileType.ELBOW_SE  # NE-corner

        # top and bottom
        for j in range(rj + 1, bj):
            if self[ri, j] == TileType.HORIZ:
                D._grid[ri, j] = TileType.BLANK
            else:  # self[ri, j] == TileType.CROSS
                D._grid[ri, j] = TileType.VERT
            if self[bi, j] == TileType.BLANK:
                D._grid[bi, j] = TileType.HORIZ
            else:  # self[bi, j] == TileType.VERT
                D._grid[bi, j] = TileType.CROSS
        # left and right
        for i in range(ri + 1, bi):
            if self[i, rj] == TileType.VERT:
                D._grid[i, rj] = TileType.BLANK
            else:  # self[i, rj] == TileType.CROSS
                D._grid[i, rj] = TileType.HORIZ
            if self[i, bj] == TileType.BLANK:
                D._grid[i, bj] = TileType.VERT
            else:  # self[i, bj] == TileType.HORIZ
                D._grid[i, bj] = TileType.CROSS
        D.rebuild()
        return D

    def huang_bump(self, a, b):
        the_cross_list = [(i, j) for (i, j) in self.all_crossings() if self.inversion_at_cross(i, j) == (a, b)]
        if len(the_cross_list) == 0:
            raise ValueError(f"No crossing found for inversion ({a}, {b})")
        x, y = the_cross_list[0]
        working_bpd = self.copy()
        working_bpd._grid[x, y] = TileType.BUMP
        return working_bpd._monk_iterate(x, y)

    def monk_insert(self, row):
        """RETURNS NORMALIZED"""
        # find easternmost tile
        working_bpd = self.normalize()
        alpha = row - 1
        col_coord = np.max(np.argwhere(self._grid[alpha, :] == TileType.ELBOW_SE))
        x, y = alpha, col_coord
        return working_bpd._monk_iterate(x, y)

    def _monk_iterate(self, x, y) -> BPD:
        working_bpd = self.resize(self.rows + 1).set_width(self.cols + 1)
        x_iter, y_iter = x, y
        while True:
            min_droops = working_bpd.min_droop_moves()
            try:
                the_move = next(((a, b), (c, d)) for ((a, b), (c, d)) in min_droops if (a == x_iter) and (b == y_iter))
            except StopIteration:
                print(f"No droop move found for position ({x_iter}, {y_iter}) in BPD:\n{working_bpd}")
                print(f"{min_droops=}")
                raise
            working_bpd = working_bpd.do_min_droop_move(the_move)
            if working_bpd[the_move[1][0], the_move[1][1]] == TileType.ELBOW_NW:
                x_iter, y_iter = the_move[1]
                continue
            z, w = the_move[1]
            # FOLLOW THE PIPE
            z_prime, w_prime = z, w
            found = False
            going_down = True
            pipe1 = working_bpd.trace_pipe(z, w, direction="down")
            pipe2 = working_bpd.trace_pipe(z, w, direction="left")
            while z_prime < working_bpd.rows and w_prime >= 0:
                tile = working_bpd[z_prime, w_prime]
                if tile == TileType.CROSS:
                    p1 = working_bpd.trace_pipe(z_prime, w_prime, direction="down")
                    p2 = working_bpd.trace_pipe(z_prime, w_prime, direction="left")
                    if (p1 == pipe1 and p2 == pipe2) or (p1 == pipe2 and p2 == pipe1):
                        found = True
                        break
                if tile == TileType.BUMP:
                    if going_down:
                        going_down = False
                        w_prime -= 1
                    else:
                        going_down = True
                        z_prime += 1
                if tile == TileType.ELBOW_SE:
                    assert not going_down
                    going_down = True
                    z_prime += 1
                if tile == TileType.ELBOW_NW:
                    assert going_down
                    going_down = False
                    w_prime -= 1
                if tile == TileType.HORIZ:
                    assert not going_down
                    w_prime -= 1
                if tile == TileType.VERT:
                    assert going_down
                    z_prime += 1
                if tile == TileType.CROSS:
                    if going_down:
                        z_prime += 1
                    else:
                        w_prime -= 1
                if tile == TileType.BUMP:
                    if going_down:
                        w_prime -= 1
                    else:
                        z_prime += 1
            if found:
                # z_prime, w_prime is a crossing
                working_bpd._grid[z, w], working_bpd._grid[z_prime, w_prime] = working_bpd._grid[z_prime, w_prime], working_bpd._grid[z, w]
                x_iter, y_iter = z_prime, w_prime
                continue
            working_bpd._grid[z, w] = TileType.CROSS
            break
        return working_bpd

    def normalize(self) -> BPD:
        if len(self) == len(self.perm):
            if len(self.perm) == self.cols:
                return self.copy()
            return self.set_width(len(self.perm))
        snap_size = len(self.perm)
        new_grid = self._grid.copy()
        if snap_size < self.cols:
            if snap_size < len(self):
                new_grid = new_grid[:snap_size, :snap_size]
            else:
                new_grid = new_grid[:, :snap_size]
        if snap_size > new_grid.shape[0] or snap_size > new_grid.shape[1]:
            new_grid = np.pad(new_grid, ((0, snap_size - new_grid.shape[0]), (0, max(0, snap_size - new_grid.shape[1]))), constant_values=TileType.TBD)
        bottom_portion = BPD.rothe_bpd(self.perm.min_coset_rep(*(list(range(self.rows)) + list(range(self.rows + 1, snap_size)))), snap_size)
        new_grid[self.rows :, :] = bottom_portion._grid[self.rows :, :]
        _invalidate_grid(new_grid)
        new_bpd = BPD(new_grid)
        assert new_bpd.rows == new_bpd.cols == snap_size
        return new_bpd

    def pop_op(self) -> tuple[BPD, tuple[int, int]]:
        # --- STEP 0 --- #
        if len(self) < len(self.perm):
            D = self.normalize()
        else:
            D = self.copy()

        # check if D has a blank tile (i.e., the coxeter length of D.w is zero)
        if self.perm.inv == 0:
            return D

        # find the first row r with a blank tile - vectorized
        blank_positions = np.argwhere(D._grid[: D.rows, : D.cols] == TileType.BLANK)
        if len(blank_positions) == 0:
            return D
        r = blank_positions[0, 0]

        # initialize the mark X at the blank tile (x,y)
        x = r
        y = np.max(blank_positions[blank_positions[:, 0] == r, 1])

        while True:
            # --- STEP 1 --- #
            # move the mark to the rightmost blank tile in the block - vectorized
            row_blanks = D._grid[x, y : D.cols] == TileType.BLANK
            if np.any(row_blanks):
                y = y + np.where(~row_blanks)[0][0] - 1 if np.any(~row_blanks) else D.cols - 1

            # --- STEP 2 --- #
            # find first j-elbow (x_,y+1) with x_>x - vectorized search
            if y + 1 >= D.cols:
                break
            elbow_positions = np.where(D._grid[x + 1 : D.rows, y + 1] == TileType.ELBOW_NW)[0]
            x_ = elbow_positions[0] + x + 1 if len(elbow_positions) > 0 else 0

            if x_ == 0:  # p==y+1
                break

            # Vectorized tile transformation for z in range(x+1, x_)
            z_range = slice(x + 1, x_)
            col_y_tiles = D._grid[z_range, y]

            # Create masks for each tile type
            blank_mask = col_y_tiles == TileType.BLANK
            horiz_mask = col_y_tiles == TileType.HORIZ
            elbow_se_mask = col_y_tiles == TileType.ELBOW_SE
            elbow_nw_mask = col_y_tiles == TileType.ELBOW_NW

            # Apply transformations
            D._grid[z_range, y][blank_mask] = TileType.VERT
            D._grid[z_range, y + 1][blank_mask] = TileType.BLANK
            D._grid[z_range, y][horiz_mask] = TileType.CROSS
            D._grid[z_range, y + 1][horiz_mask] = TileType.HORIZ
            D._grid[z_range, y][elbow_se_mask] = TileType.VERT
            D._grid[z_range, y + 1][elbow_se_mask] = TileType.ELBOW_SE
            D._grid[z_range, y][elbow_nw_mask] = TileType.CROSS
            D._grid[z_range, y + 1][elbow_nw_mask] = TileType.ELBOW_NW

            D._grid[x, y] = TileType.ELBOW_SE  # NW-corner
            D._grid[x_, y + 1] = TileType.BLANK  # SE-corner
            if D[x_, y] == TileType.ELBOW_SE:  # SW-corner
                D._grid[x_, y] = TileType.VERT
            else:  # D._grid[x_,y] == TileType.HORIZ
                D._grid[x_, y] = TileType.ELBOW_NW
            if D[x, y + 1] == TileType.ELBOW_SE:  # NE-corner
                D._grid[x, y + 1] = TileType.HORIZ
            else:  # D._grid[x,y+1] == TileType.VERT
                D._grid[x, y + 1] = TileType.ELBOW_NW

            # move the mark X to the SE-corner of U
            x = x_
            y = y + 1

        # --- STEP 3 --- #
        a = y  # where (x,y) is the final position of the mark X

        # Find x_ - vectorized search from bottom
        if y < D.cols:
            elbow_positions = np.where(D._grid[x + 1 : D.rows, y] == TileType.ELBOW_SE)[0]
            x_ = elbow_positions[-1] + x + 1 if len(elbow_positions) > 0 else 0
        else:
            x_ = 0

        # Vectorized tile transformation for z in range(x+1, x_)
        if x_ > x + 1:
            z_range = slice(x + 1, x_)
            col_y_tiles = D._grid[z_range, y]

            # Create masks for each tile type
            blank_mask = col_y_tiles == TileType.BLANK
            horiz_mask = col_y_tiles == TileType.HORIZ
            elbow_se_mask = col_y_tiles == TileType.ELBOW_SE
            elbow_nw_mask = col_y_tiles == TileType.ELBOW_NW

            # Apply transformations
            D._grid[z_range, y][blank_mask] = TileType.VERT
            D._grid[z_range, y + 1][blank_mask] = TileType.BLANK
            D._grid[z_range, y][horiz_mask] = TileType.CROSS
            D._grid[z_range, y + 1][horiz_mask] = TileType.HORIZ
            D._grid[z_range, y][elbow_se_mask] = TileType.VERT
            D._grid[z_range, y + 1][elbow_se_mask] = TileType.ELBOW_SE
            D._grid[z_range, y][elbow_nw_mask] = TileType.CROSS
            D._grid[z_range, y + 1][elbow_nw_mask] = TileType.ELBOW_NW

        if x_ > 0:
            D._grid[x, y] = TileType.ELBOW_SE  # NW-corner
            D._grid[x_, y + 1] = TileType.ELBOW_SE  # SE-corner
            D._grid[x_, y] = TileType.VERT  # SW-corner
            if D[x, y + 1] == TileType.ELBOW_SE:  # NE-corner
                D._grid[x, y + 1] = TileType.HORIZ
            else:  # D._grid[x,y+1] == TileType.VERT
                D._grid[x, y + 1] = TileType.ELBOW_NW

        return D.resize(self.rows), (a + 1, r + 1)

    def column_perm_at_row(self, row: int) -> Permutation:
        build_perm = []
        for col in range(self.cols):
            if self[row, col].entrance_from_bottom:
                try:
                    build_perm.append(self.trace_pipe(row, col))
                except ValueError:
                    build_perm.append(self.trace_pipe(row, col, direction="down"))
            else:
                build_perm.append(None)
        return Permutation.from_partial(build_perm)

    def resize(self, new_num_rows: int, check_valid=True) -> BPD:
        if new_num_rows > self.rows:
            new_bpd = self.copy()
            if check_valid:
                new_bpd = self.normalize()
            if new_bpd.rows < new_num_rows:
                new_bpd._grid = np.pad(new_bpd._grid, ((0, new_num_rows - new_bpd.rows), (0, max(0, new_num_rows - new_bpd.cols))), constant_values=TileType.TBD)
            elif new_bpd.rows > new_num_rows:
                new_bpd._grid = new_bpd._grid[:new_num_rows, : max(new_num_rows, len(self.perm))]
            new_bpd.rebuild()
            if check_valid:
                assert new_bpd.is_valid, f"Resulting BPD is not valid after increasing size, \n{pretty(self)} {new_num_rows} {new_bpd!r}"
            return new_bpd
        if check_valid and new_num_rows < len(self.perm):
            new_bpd = BPD(self._grid[:new_num_rows, :])
            new_bpd.rebuild()
            return new_bpd
        return BPD(self._grid[:new_num_rows, :])

    @classmethod
    def from_rc_graph(cls, rc_graph) -> BPD:
        if cls.DEBUG:
            assert len(rc_graph.perm.trimcode) <= len(rc_graph), f"RC graph permutation length exceeds RC graph size: {rc_graph.perm} vs {len(rc_graph)}"
            assert rc_graph.perm.inv == sum(len(rc_graph[i]) for i in range(len(rc_graph))), f"RC graph permutation inversion count does not match RC graph crossings: {rc_graph.perm} vs {rc_graph}"
        num_rows = len(rc_graph)
        n = max(num_rows, len(rc_graph.perm))
        bpd = BPD(np.full((n, n), fill_value=TileType.TBD, dtype=TileType))
        coords = [rc_graph.left_to_right_inversion_coords(i) for i in range(rc_graph.perm.inv - 1, -1, -1)]

        bpd = bpd.inverse_pop_op(*[(i + j - 1, i) for i, j in coords])

        if cls.DEBUG:
            assert bpd.perm.inv == len(bpd.all_blanks())
            assert bpd.perm == rc_graph.perm
        if len(bpd) != len(rc_graph):
            return bpd.resize(num_rows)
        return bpd

    def combine(self, other, shift=None) -> BPD:
        """Shift the BPD up by adding empty rows at the bottom."""
        return BPD.from_rc_graph(RCGraph([*self.to_rc_graph()[:shift], *other.to_rc_graph().shiftup(shift)]))

    # def complete_partial(self, partition = None) -> BPD:
    #     """Complete a partial BPD to a full BPD."""
    #     if partition is None:
    #         partition = [self.cols] * self.rows
    #     if len(partition) < len(self.perm):
    #         partition = partition + [0] * (len(self.perm) - len(partition))
    #     new_grid = self._grid.copy()
    #     new_grid = np.pad(new_grid, ((0, self.rows - len(partition)), (0, max(0, len(self.perm) - self.cols))), constant_values=TileType.TBD)

    # def left_row_act(self, p: int):
    #     new_grid = self._grid.copy()
    #     new_grid = np.pad(new_grid, ((1, 0), (0,1)), constant_values=TileType.TBD)

    def product(self, other: BPD) -> dict[BPD, int]:
        """Compute the product of this BPD with another."""
        from schubmult.utils.perm_utils import add_perm_dict

        if len(other) == 0:
            return {self: 1}
        self_perm = self.perm
        other_shift = other.shiftup(len(self))
        if self_perm.inv == 0:
            # return {BPD.rothe_bpd(Permutation([]), len(self) + len(other)).inverse_pop_op(*other_reduced_compatible).resize(len(self) + len(other)): 1}
            return {other_shift: 1}
        self_len = len(self)
        other_perm_len = len(other.perm)
        num_zeros = max(len(other), other_perm_len)
        if self.DEBUG:
            assert len(self_perm.trimcode) <= self_len, f"{self=}, {self_perm=}"
        base_bpd = self
        # print("BASE")
        # print(base_bpd)
        buildup_module = {base_bpd: 1}

        for _ in range(num_zeros):
            new_buildup_module = {}
            # print("FATPANT")
            for bpd, coeff in buildup_module.items():
                # print("BOIP")
                # print(bpd)
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(bpd.right_zero_act(), coeff))
            buildup_module = new_buildup_module
        ret_module = {}
        other_len = other.rows
        self_perm_inv = self_perm.inv
        other_perm_inv = other.perm.inv
        target_inv = self_perm_inv + other_perm_inv
        for bpd, coeff in buildup_module.items():
            # print(bpd)

            if self.DEBUG:
                assert len(bpd.perm.trimcode) <= len(bpd), f"{bpd=}, {bpd.perm=}"
            if (bpd.perm * other_shift.perm).inv != bpd.perm.inv + other_perm_inv:
                continue
            new_bpd = other_shift.inverse_pop_op(*bpd.as_reduced_compatible())
            new_bpd_perm = new_bpd.perm
            if len(new_bpd_perm.trimcode) > other_len + self_len:
                continue
            new_bpd = new_bpd.resize(other_len + self_len).snap_width()
            if new_bpd_perm != new_bpd.perm:
                continue
            if new_bpd.is_valid and new_bpd.is_reduced and new_bpd_perm.inv == target_inv and len(new_bpd) <= self_len + other_len:
                ret_module = add_perm_dict(ret_module, {new_bpd: coeff})

        return ret_module

    def prod_with_bpd(self, other: BPD) -> BPD:
        """Deprecated: Use product() instead. Returns the single BPD from product dictionary."""
        return self.product(other)

    def inverse_pop_op(self, *interlaced_rc) -> BPD:
        if self.rows < len(self.perm):
            D = self.normalize()
        else:
            D = self.copy()

        if len(interlaced_rc) == 2 and isinstance(interlaced_rc[0], int):
            interlaced_rc = [interlaced_rc]
        else:
            interlaced_rc = list(reversed(interlaced_rc))

        while len(interlaced_rc) > 0:
            a, r = interlaced_rc.pop()
            if D.rows <= a or D.rows <= r:
                new_num_rows = max(D.rows, a + 1, r + 1)
                D = D.resize(new_num_rows)
            # find first elbow in column a - vectorized search
            elbow_positions = np.where(D._grid[: D.rows, a] == TileType.ELBOW_SE)[0]
            if len(elbow_positions) == 0:
                raise ValueError("No elbow found in specified column for inverse pop operation when inserting ")
            x_ = elbow_positions[-1]  # Last (bottom-most) elbow
            D._grid[x_, a] = TileType.CROSS
            y = a - 1
            # find x in column y, it will be an SE elbow - vectorized search
            if y >= 0:
                elbow_positions = np.where(D._grid[:x_, y] == TileType.ELBOW_SE)[0]
                if len(elbow_positions) == 0:
                    raise ValueError("No elbow found in specified column for inverse pop operation")
                x = elbow_positions[-1]  # Last (bottom-most) elbow before x_
            else:
                x = -1
                raise ValueError("No elbow found in specified column for inverse pop operation")
            D._grid[x, y] = TileType.BLANK

            # Vectorized tile transformation for z in range(x+1, x_)
            if x_ > x + 1:
                z_range = slice(x + 1, x_)
                col_y = D._grid[z_range, y]
                col_y1 = D._grid[z_range, y + 1]

                # Create masks for tile patterns
                mask1 = (col_y == TileType.VERT) & (col_y1 == TileType.BLANK)
                mask2 = (col_y == TileType.CROSS) & (col_y1 == TileType.ELBOW_NW)
                mask3 = (col_y == TileType.CROSS) & (col_y1 == TileType.HORIZ)
                mask4 = (col_y == TileType.VERT) & (col_y1 == TileType.ELBOW_SE)

                # Apply transformations
                D._grid[z_range, y][mask1] = TileType.BLANK
                D._grid[z_range, y + 1][mask1] = TileType.VERT
                D._grid[z_range, y][mask2] = TileType.TBD
                D._grid[z_range, y + 1][mask2] = TileType.TBD
                D._grid[z_range, y][mask3] = TileType.TBD
                D._grid[z_range, y + 1][mask3] = TileType.CROSS
                D._grid[z_range, y][mask4] = TileType.ELBOW_SE
                D._grid[z_range, y + 1][mask4] = TileType.CROSS

            D.rebuild()

            while True:
                # --- STEP 1 --- #
                # move the mark to the rightmost blank tile in the block
                if x == r - 1:
                    break

                # find x_
                x_ = x
                y = y - 1
                if self.DEBUG:
                    assert D[x_, y + 1] == TileType.BLANK, "Expected NW elbow during inverse pop operation"
                # Vectorized search: find first non-blank from left in row
                if y >= 0:
                    row_non_blanks = D._grid[x_, : y + 1] != TileType.BLANK
                    non_blank_positions = np.where(row_non_blanks)[0]
                    if len(non_blank_positions) > 0:
                        y = non_blank_positions[-1]  # Rightmost non-blank
                    else:
                        y = -1
                else:
                    y = -1
                D._grid[x_, y + 1] = TileType.ELBOW_NW

                # find x at SE elbow - vectorized search
                if y >= 0:
                    elbow_positions = np.where(D._grid[:x_, y] == TileType.ELBOW_SE)[0]
                    if len(elbow_positions) == 0:
                        raise ValueError("No elbow found in specified column for inverse pop operation")
                    x = elbow_positions[-1]
                else:
                    raise ValueError("No elbow found in specified column for inverse pop operation")

                # [x, y] becomes
                D._grid[x, y] = TileType.BLANK

                # Vectorized tile transformation for z in range(x+1, x_)
                if x_ > x + 1:
                    z_range = slice(x + 1, x_)
                    col_y = D._grid[z_range, y]
                    col_y1 = D._grid[z_range, y + 1]

                    # Create masks for tile patterns
                    mask1 = (col_y == TileType.VERT) & (col_y1 == TileType.BLANK)
                    mask2 = (col_y == TileType.CROSS) & (col_y1 == TileType.ELBOW_NW)
                    mask3 = (col_y == TileType.CROSS) & (col_y1 == TileType.HORIZ)
                    mask4 = (col_y == TileType.VERT) & (col_y1 == TileType.ELBOW_SE)

                    # Apply transformations
                    D._grid[z_range, y][mask1] = TileType.BLANK
                    D._grid[z_range, y + 1][mask1] = TileType.VERT
                    D._grid[z_range, y][mask2] = TileType.TBD
                    D._grid[z_range, y + 1][mask2] = TileType.TBD
                    D._grid[z_range, y][mask3] = TileType.TBD
                    D._grid[z_range, y + 1][mask3] = TileType.CROSS
                    D._grid[z_range, y][mask4] = TileType.ELBOW_SE
                    D._grid[z_range, y + 1][mask4] = TileType.CROSS

                D.rebuild()
        return D

    @property
    def is_reduced(self):
        if not self.is_valid:
            return False
        return self.perm.inv == np.sum(self.length_vector)

    @cache
    def as_reduced_compatible(self):
        if not self.is_valid:
            raise ValueError("BPD is not valid, cannot compute reduced compatible sequence")
        if not self.is_reduced:
            raise ValueError("BPD is not reduced, cannot compute reduced compatible sequence")
        work_bpd = self
        ret = []
        last_inv = work_bpd.perm.inv
        while work_bpd.perm.inv > 0:
            work_bpd, (reflection, row) = work_bpd.pop_op()
            if self.DEBUG:
                assert work_bpd.perm.inv == last_inv - 1, "Invariant violated in as_reduced_compatible"
            last_inv = work_bpd.perm.inv
            ret.append((reflection, row))
        ret.reverse()
        return tuple(ret)

    def rebuild(self) -> None:
        """Rebuild the BPD to resolve any TBD tiles"""
        # Keep only BLANK and CROSS tiles, wipe out everything else to TBD
        _invalidate_grid(self._grid)
        self._perm = None
        self._valid = None
        self._word = None
        self._unzero_cache = None
        self.build()

    def zero_out_last_row(self) -> BPD:
        return self.resize(self.rows - 1)

    def set_tile(self, i: int, j: int, tile_type: TileType) -> None:
        new_bpd = self.copy()
        new_bpd._grid[i, j] = tile_type
        new_bpd.rebuild()
        return new_bpd

    def right_zero_act(self) -> set[BPD]:
        # # find crosses, untransition them
        if self._unzero_cache is not None:
            return self._unzero_cache
        if not self.is_valid:
            return set()
        resized = self.resize(self.rows + 1)
        results = set()
        min_width = len(self.perm) + 1
        if resized.cols < min_width:
            resized = resized.set_width(min_width)
        results = set()

        crossings = np.array([(i, j) for (i, j) in resized.all_crossings() if i == self.rows], dtype=np.int32)

        n_crossings = len(crossings)
        base_grid = resized._grid.copy()
        for mask in range(1 << n_crossings):  # 2^n combinations
            # Start from base grid instead of full BPD copy
            working_grid = base_grid.copy()
            # Apply mask: set to TBD where mask bit is 1
            for idx in range(n_crossings):
                if mask & (1 << idx):
                    i, j = crossings[idx]
                    working_grid[i, j] = TileType.TBD
            _invalidate_grid(working_grid)
            new_bpd = BPD(working_grid).snap_width()

            if new_bpd.is_valid and new_bpd.is_reduced:
                results.add(new_bpd)
        results = set()

        # up_perms = [perm for perm, _ in ASx(self.perm, self.rows) * ASx([], 1)]

        # self_perm = self.perm
        # self_snap = self.snap_width()
        # for up_perm in up_perms:
        #     # want
        #     if up_perm == self_perm:
        #         results.add(self_snap.resize(self.rows + 1))
        #         continue
        #     new_grid = self_snap._grid.copy()
        #     row_pad = 1
        #     if len(up_perm) > len(self_perm):
        #         col_pad = len(up_perm) - len(self_perm)
        #     new_grid = np.pad(new_grid, ((0, row_pad), (0, col_pad)), constant_values=TileType.TBD)
        #     moving_perm =

        self._unzero_cache = results
        return results

    def set_tiles(self, a, b, value: TileType) -> None:
        ret = self.copy()
        ret._grid[a, b] = value
        ret.rebuild()
        return ret

    #        for p in range(len(v) + 1 - vnum):
    #     vpm_list2 = []
    #     for vpm, b in vpm_list:
    #         if vpm[vnum - 1] == len(v) + 1:
    #             vpm2 = [*vpm]
    #             vpm2.pop(vnum - 1)
    #             vp = permtrim(vpm2)
    #             ret_list.add(
    #                 (
    #                     tuple([v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))]),
    #                     vp,
    #                 ),
    #             )
    #         for j in range(vnum, len(vup) + 2):
    #             if vpm[j] <= b:
    #                 continue
    #             for i in range(vnum):
    #                 if has_bruhat_ascent(vpm, i, j):
    #                     vpm_list2 += [(vpm.swap(i, j), vpm[j])]
    #     vpm_list = vpm_list2
    # for vpm, b in vpm_list:
    #     if vpm[vnum - 1] == len(v) + 1:
    #         vpm2 = [*vpm]
    #         vpm2.pop(vnum - 1)
    #         vp = permtrim(vpm2)
    #         ret_list.add(
    #             (
    #                 tuple([v[i] for i in range(vnum, len(v)) if ((i > len(vp) and v[i] == i) or (i <= len(vp) and v[i] == vp[i - 1]))]),
    #                 vp,
    #             ),
    #         )

    # n_crossings = len(crossings)
    # base_grid = resized._grid.copy()
    # for mask in range(1 << n_crossings):  # 2^n combinations
    #     # Start from base grid instead of full BPD copy
    #     working_grid = base_grid.copy()
    #     # Apply mask: set to TBD where mask bit is 1
    #     for idx in range(n_crossings):
    #         if mask & (1 << idx):
    #             i, j = crossings[idx]
    #             working_grid[i, j] = TileType.TBD
    #     _invalidate_grid(working_grid)
    #     new_bpd = BPD(working_grid).snap_width()

    #     if new_bpd.is_valid and new_bpd.is_reduced:
    #         results.add(new_bpd)

    # return results

    def snap_width(self) -> BPD:
        """Snap the width of the BPD to the length of its permutation."""
        perm_length = len(self.perm)
        if self.cols == perm_length or (self.rows > perm_length and self.cols == self.rows):
            return self
        return self.set_width(max(self.rows, perm_length))

    def polyvalue(self, x: Sequence[Expr], y: Sequence[Expr] | None = None, **_kwargs) -> Expr:
        """
        Compute the Schubert polynomial value for this BPD.

        Args:
            x: Variable or list of variables for polynomial
            y: Optional second set of variables for double Schubert polynomial
            **_kwargs: Additional keyword arguments for polynomial computation (unused)
        """
        from schubmult.symbolic import prod

        if y is None:
            return prod(x[i + 1] for i, _ in self.all_blanks())
        return prod(x[i + 1] - y[j + 1] for i, j in self.all_blanks())

    def to_rc_graph(self) -> RCGraph:
        """
        Convert this BPD to an RC-graph representation.

        Returns:
            RCGraph object (if available in the module)
        """

        rows = [[] for _ in range(self.rows)]
        rc = self.as_reduced_compatible()
        for reflection, row in reversed(rc):
            rows[row - 1] = rows[row - 1] + [reflection]

        return RCGraph([tuple(r) for r in rows])

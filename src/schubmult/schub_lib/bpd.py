"""
Bumpless Pipe Dreams (BPD) module

A bumpless pipe dream is a pipe dream diagram with no elbows - only crossings and empty boxes.
Represented as an n x n grid where:
- 1 = crossing (pipes cross)
- 0 = empty box (pipes go straight through)

For general pipe dreams, there are 6 possible tile types.
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
        return self in (TileType.HORIZ, TileType.ELBOW_SE, TileType.CROSS)

    @cached_property
    def feeds_up(self) -> bool:
        """True if the vertical pipe continues upwards"""
        return self in (TileType.VERT, TileType.ELBOW_NW, TileType.CROSS)

    @cached_property
    def entrance_from_bottom(self) -> bool:
        """True if a pipe can enter from the bottom"""
        return self in (TileType.VERT, TileType.ELBOW_SE, TileType.CROSS)

    @cached_property
    def entrance_from_left(self) -> bool:
        """True if a pipe can enter from the left"""
        return self in (TileType.HORIZ, TileType.ELBOW_NW, TileType.CROSS)


def _invalidate_grid(grid: np.ndarray) -> None:
    """Helper function to invalidate a grid by setting TBD tiles"""
    # Compute mask once: keep only BLANK and CROSS tiles
    mask = (grid == TileType.BLANK) | (grid == TileType.CROSS)
    grid[~mask] = TileType.TBD


def _get_northeast_pipe_count(grid, target_row, target_col):
    """
    Computes the number of pipes weakly northeast of a crossing at (target_row, target_col).
    Grid is a 2D list of TileType.
    """
    # r[i][j] stores the corner sum r_A(i, j)
    # Using n+1 to handle 0-th row/col boundary conditions easily
    r = np.zeros((grid.shape[0] + 1, grid.shape[1] + 1), dtype=int)

    for i in range(1, grid.shape[0] + 1):
        for j in range(1, grid.shape[1] + 1):
            tile = grid[i - 1, j - 1]

            # Based on Lemma 3.5 and 3.8:
            # r(i,j) - r(i-1, j-1) is:
            # 0 if BLANK, 2 if CROSS, 1 otherwise.
            nw_val = r[i - 1, j - 1]

            if tile == TileType.BLANK:
                diff = 0
            elif tile == TileType.CROSS:
                diff = 2
            else:
                diff = 1

            r[i, j] = nw_val + diff

    # The label is r_A(i, j) - 1.
    # Therefore, the total count weakly northeast is simply r_A(i, j).
    return r[target_row + 1, target_col + 1]


def _is_asm(arr):
    # Ensure it's a square numpy array
    n, m = arr.shape
    if n != m:
        return False

    # 1. Check Row and Column Sums (Must be exactly 1)
    if not np.all(arr.sum(axis=0) == 1) or not np.all(arr.sum(axis=1) == 1):
        return False

    # 2. Check Allowed Values (Only 0, 1, -1)
    # Using absolute value is a fast way to check this
    # if not np.all(np.abs(arr) <= 1):
    #     return False

    # 3. Check Alternating Signs and Boundary +1s
    # In a valid ASM, the cumulative sum along any row/column
    # must stay between 0 and 1 at all times.
    # Because row sums are 1, the first/last non-zero must be 1
    # if the prefix sums are always 0 or 1.

    # Check rows
    row_cumsum = np.cumsum(arr, axis=1)
    if np.any((row_cumsum < 0) | (row_cumsum > 1)):
        return False

    # Check columns
    col_cumsum = np.cumsum(arr, axis=0)
    if np.any((col_cumsum < 0) | (col_cumsum > 1)):
        return False

    # 4. Final verification: The non-zero entries must actually change.
    # Example of a failure mode for cumsum alone: [1, 0, 0] is fine,
    # but [1, 1, -1] has cumsums [1, 2, 1], which is caught above.
    # The only thing left is to ensure we don't have two 1s or two -1s in a row.
    # We can do this by checking if the absolute difference of non-zero entries is 2.

    for axis in [0, 1]:
        # Flattened check for each row/column
        for line in arr if axis == 1 else arr.T:
            non_zeros = line[line != 0]
            # Since we know they sum to 1 and cumsum is 0/1,
            # we just need to ensure no two adjacent non-zeros are the same.
            if np.any(np.diff(non_zeros) == 0):
                return False

    return True


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

        # Create lookup table: shape (9, 9) for tile types 0-7 plus None (-1 -> index 8)
        lookup = np.full((9, 9), TileType.ELBOW_SE, dtype=TileType)

        # Map None to index 8
        none_idx = 8

        # For each combination of (left_tile, up_tile), compute result
        for left_val in range(-1, 8):  # -1 for None, 0-7 for tile types
            left_idx = none_idx if left_val == -1 else left_val
            left_feeds_right = False if left_val == -1 else TileType(left_val).feeds_right if left_val < 8 else False

            for up_val in range(-1, 8):
                up_idx = none_idx if up_val == -1 else up_val
                up_entrance = False if up_val == -1 else TileType(up_val).entrance_from_bottom if up_val < 8 else False

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

        good_cols = []
        for col in range(self.cols):
            if self[self.rows - 1, col].entrance_from_bottom:
                good_cols.append(col + 1)
            elif self[self.rows - 1, col] == TileType.TBD:
                raise ValueError("Cannot compute permutation with unresolved TBD tiles")
        # if len(good_cols) == 0:
        #     self._perm = Permutation([])
        #     print(self)
        #     return self._perm

        small_perm = Permutation.ref_product(*self.word)
        build_perm = [good_cols[small_perm[i] - 1] if small_perm[i] - 1 < len(good_cols) else small_perm[i] - 1 for i in range(len(good_cols))] + [None] * (
            max(good_cols, default=0) + 1 - len(good_cols)
        )
        self._perm = self._column_perm * Permutation.from_partial(build_perm)
        return self._perm
        # if self._perm is not None:
        #     return self._perm
        # buildperm = []

        # # Cache property access to avoid overhead in tight loop
        # rows = self._grid.shape[0]
        # cols = self._grid.shape[1]
        # grid = self._grid  # Direct array access

        # for row in range(rows):
        #     # trace from the right
        #     going_left = True
        #     current_row = row
        #     current_col = cols - 1
        #     while current_col >= 0 and current_row < rows:
        #         tile = grid[current_row, current_col]
        #         if tile == TileType.ELBOW_NW:
        #             # Elbow NW: go down
        #             assert not going_left, "Invalid pipe direction at ELBOW_NW"
        #             current_col -= 1
        #             going_left = True
        #         elif tile == TileType.ELBOW_SE:
        #             # Elbow SE: go right
        #             assert going_left, "Invalid pipe direction at ELBOW_SE"
        #             current_row += 1
        #             going_left = False
        #         elif tile == TileType.CROSS:
        #             if going_left:
        #                 current_col -= 1
        #             else:
        #                 current_row += 1
        #         elif tile == TileType.HORIZ:
        #             assert going_left
        #             current_col -= 1
        #         elif tile == TileType.VERT:
        #             assert not going_left
        #             current_row += 1
        #         else:
        #             raise ValueError(f"Invalid tile type {tile} in BPD")
        #     buildperm.append(current_col + 1)
        # n = max(buildperm, default=1)
        # self._perm = self._column_perm * Permutation.from_partial(buildperm + [None] * (n - len(buildperm)))
        # return self._perm

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
        # asm = np.zeros(shape=(len(perm), len(perm)), dtype=np.int8)
        # graph = [(i - 1, j - 1) for (i, j) in perm.graph]
        # rows, cols = zip(*graph) if graph else ([], [])
        # asm[rows, cols] = 1
        # bpd = cls.from_asm(asm)
        # # print(asm)
        # # print(bpd)
        # if num_rows is None or num_rows == len(perm):
        #     return bpd
        # return bpd.resize(num_rows)

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

        # Read crossings from bottom to top, left to right
        for col in range(self.cols):
            for row in range(self.rows - 1, -1, -1):
                if self[row, col] == TileType.CROSS:
                    # Count pipes weakly northeast of this crossing
                    # Weakly northeast: r <= row and c >= col

                    word.append(_get_northeast_pipe_count(self._grid, row, col) - 1)

        self._word = tuple(word)
        return self._word

    @property
    def reduced_word(self) -> Tuple[int, ...]:
        """
         Compute the reduced word by reading crossings row by row from top to bottom.

        For each crossing at position (i,j), the word value is the number of pipes
         weakly northeast of the crossing minus 1.

         Returns:
             Tuple of integers representing a reduced word
        """
        word = []

        # Read from top to bottom, left to right
        for row in range(self.rows):
            for col in range(self.cols):
                if self[row, col] == 1:
                    # Count pipes weakly northeast of this crossing
                    pipes_northeast = self.rows - col
                    word.append(pipes_northeast)

        return tuple(word)

    def set_width(self, width):
        """Set the width of the BPD by adding empty columns on the right if needed."""
        if width < self.cols:
            if len(self.perm) > width:
                raise ValueError("New width must be at least the length of the permutation")
            new_grid = np.full((self.rows, width), fill_value=TileType.TBD, dtype=TileType)
            new_grid[:, :width] = self._grid[:, :width]
            bop = BPD(new_grid)
            bop.rebuild()
            assert bop.is_valid, f"Resulting BPD is not valid after reducing width, {pretty(self)} {width} "
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
        self._valid = _is_asm(self.to_asm())
        return self._valid

        # Cache shape to avoid property access overhead
        # rows = self._grid.shape[0]
        # cols = self._grid.shape[1]

        # # Left column boundary check - vectorized
        # left_col = self._grid[:rows, 0]
        # if np.any((left_col == TileType.ELBOW_NW) | (left_col == TileType.CROSS) | (left_col == TileType.HORIZ)):
        #     return False

        # # Right column boundary check - vectorized
        # right_col = self._grid[:rows, cols - 1]
        # if np.any((right_col == TileType.ELBOW_NW) | (right_col == TileType.VERT)):
        #     return False

        # # Top row boundary check - vectorized
        # top_row = self._grid[0, :cols]
        # if np.any((top_row == TileType.ELBOW_NW) | (top_row == TileType.VERT) | (top_row == TileType.CROSS)):
        #     return False

        # # Interior connectivity checks - vectorized using boolean masks
        # if rows > 2 and cols > 2:
        #     # Get interior cells
        #     interior = self._grid[1 : rows - 1, 1 : cols - 1]
        #     right_neighbors = self._grid[1 : rows - 1, 2:cols]
        #     left_neighbors = self._grid[1 : rows - 1, 0 : cols - 2]
        #     up_neighbors = self._grid[0 : rows - 2, 1 : cols - 1]
        #     down_neighbors = self._grid[2:rows, 1 : cols - 1]

        #     # Check right connectivity: feeds_right requires entrance_from_left
        #     feeds_right = (interior == TileType.HORIZ) | (interior == TileType.ELBOW_SE) | (interior == TileType.CROSS)
        #     entrance_left = (right_neighbors == TileType.HORIZ) | (right_neighbors == TileType.ELBOW_NW) | (right_neighbors == TileType.CROSS)
        #     if np.any(feeds_right & ~entrance_left):
        #         return False

        #     # Check up connectivity: feeds_up requires entrance_from_bottom
        #     feeds_up = (interior == TileType.VERT) | (interior == TileType.ELBOW_NW) | (interior == TileType.CROSS)
        #     entrance_bottom = (up_neighbors == TileType.VERT) | (up_neighbors == TileType.ELBOW_SE) | (up_neighbors == TileType.CROSS)
        #     if np.any(feeds_up & ~entrance_bottom):
        #         return False

        #     # Check left connectivity: entrance_from_left requires feeds_right
        #     entrance_left_interior = (interior == TileType.HORIZ) | (interior == TileType.ELBOW_NW) | (interior == TileType.CROSS)
        #     feeds_right_left = (left_neighbors == TileType.HORIZ) | (left_neighbors == TileType.ELBOW_SE) | (left_neighbors == TileType.CROSS)
        #     if np.any(entrance_left_interior & ~feeds_right_left):
        #         return False

        #     # Check bottom connectivity: entrance_from_bottom requires feeds_up
        #     entrance_bottom_interior = (interior == TileType.VERT) | (interior == TileType.ELBOW_SE) | (interior == TileType.CROSS)
        #     feeds_up_down = (down_neighbors == TileType.VERT) | (down_neighbors == TileType.ELBOW_NW) | (down_neighbors == TileType.CROSS)
        #     if np.any(entrance_bottom_interior & ~feeds_up_down):
        #         return False

        # try:
        #     return self.perm.inv == sum(self.length_vector)
        # except Exception:
        #     return False

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
        return set(zip(*np.where(self._grid == tile_type)))

    def droop_moves(self) -> set[tuple[tuple[int, int], tuple[int, int]]]:
        import itertools

        droop_moves = set()
        for (ri, rj), (bi, bj) in itertools.product(self.all_se_elbows(), self.all_blanks()):
            if bi > ri and bj > rj:
                LEGAL = True
                for i, j in itertools.product(range(ri, bi + 1), range(rj, bj + 1)):
                    if (i, j) != (ri, rj) and (self[i, j] == TileType.ELBOW_SE or self[i, j] == TileType.ELBOW_NW):  # if not NW-corner, check if elbow
                        LEGAL = False
                        break
                if LEGAL:
                    droop_moves.add(((ri, rj), (bi, bj)))
        return droop_moves

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

    def resize(self, new_num_rows: int, column_perm: Permutation = Permutation([])) -> BPD:
        if new_num_rows > self.rows:
            new_bpd = self.normalize()
            if new_bpd.rows < new_num_rows:
                new_bpd._grid = np.pad(new_bpd._grid, ((0, new_num_rows - new_bpd.rows), (0, max(0, new_num_rows - self.cols))), constant_values=TileType.TBD)
            elif new_bpd.rows > new_num_rows:
                new_bpd._grid = new_bpd._grid[:new_num_rows, : max(new_num_rows, len(self.perm))]
            new_bpd.rebuild()
            assert new_bpd.is_valid, f"Resulting BPD is not valid after increasing size, \n{pretty(self)} {new_num_rows} "
            return new_bpd
        if new_num_rows < len(self.perm):
            return BPD(self._grid[:new_num_rows, :], column_perm=self.column_perm_at_row(new_num_rows - 1) if column_perm is None else column_perm)
        return BPD(self._grid[:new_num_rows, :], column_perm=column_perm)

    @classmethod
    def from_rc_graph(cls, rc_graph) -> BPD:
        assert len(rc_graph.perm.trimcode) <= len(rc_graph), f"RC graph permutation length exceeds RC graph size: {rc_graph.perm} vs {len(rc_graph)}"
        num_rows = len(rc_graph)
        n = max(num_rows, len(rc_graph.perm))
        bpd = BPD(np.full((n, n), fill_value=TileType.TBD, dtype=TileType))
        coords = [rc_graph.left_to_right_inversion_coords(i) for i in range(rc_graph.perm.inv - 1, -1, -1)]

        bpd = bpd.inverse_pop_op(*[(i + j - 1, i) for i, j in coords])

        assert bpd.perm.inv == len(bpd.all_blanks())
        assert bpd.perm == rc_graph.perm
        if len(bpd) != len(rc_graph):
            return bpd.resize(num_rows)
        return bpd

    def product(self, other: BPD) -> dict[BPD, int]:
        """Compute the product of this BPD with another."""
        from schubmult.utils.perm_utils import add_perm_dict

        other_graph = other.to_rc_graph()
        assert len(other.perm.trimcode) <= len(other_graph), f"{other=}, {other_graph=}"
        # other_reduced_compatible = [(a + len(self), r + len(self)) for a, r in other.as_reduced_compatible()]
        # other_reduced_compatible.reverse()
        self_perm = self.perm
        if self_perm.inv == 0:
            # return {BPD.rothe_bpd(Permutation([]), len(self) + len(other)).inverse_pop_op(*other_reduced_compatible).resize(len(self) + len(other)): 1}
            return {BPD.from_rc_graph(other_graph.prepend(len(self))): 1}
        self_len = len(self)
        other_perm_len = len(other.perm)
        num_zeros = max(len(other), other_perm_len)
        assert len(self_perm.trimcode) <= self_len, f"{self=}, {self_perm=}"
        base_bpd = self
        buildup_module = {base_bpd: 1}

        for _ in range(num_zeros):
            new_buildup_module = {}
            for bpd, coeff in buildup_module.items():
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(bpd.right_zero_act(), coeff))
            buildup_module = new_buildup_module
        ret_module = {}

        self_perm_inv = self_perm.inv
        other_perm_inv = other.perm.inv
        target_inv = self_perm_inv + other_perm_inv
        for bpd, coeff in buildup_module.items():
            assert bpd.is_valid, f"Invalid BPD in product buildup: {pretty(bpd)}"
            new_rc = RCGraph([*bpd.to_rc_graph()[:self_len], *other_graph.shiftup(self_len)])
            if new_rc.is_valid and len(new_rc.perm.trimcode) <= len(new_rc):
                new_bpd = BPD.from_rc_graph(new_rc)
                assert len(new_bpd) == self_len + len(other)

                new_bpd_perm = new_bpd.perm
                if new_bpd.is_valid and new_bpd_perm.inv == target_inv and len(new_bpd_perm.trimcode) <= len(new_bpd):
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
                # replace with NW elbow
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
        return self.perm.inv == sum(self.length_vector)

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
        self.build()
        return self.is_valid

    def zero_out_last_row(self) -> BPD:
        return self.resize(self.rows - 1, column_perm=Permutation([]))

    def set_tile(self, i: int, j: int, tile_type: TileType) -> None:
        new_bpd = self.copy()
        new_bpd._grid[i, j] = tile_type
        new_bpd.rebuild()
        return new_bpd

    # def __setitem__(self, key: tuple[int, int], value: TileType) -> None:
    #     i, j = key
    #     self._grid[i, j] = value
    #     self.rebuild()

    def right_zero_act(self) -> set[BPD]:
        # # find crosses, untransition them
        # from schubmult import ASx
        if not self.is_valid:
            return set()
        resized = self.resize(self.rows + 1)
        min_width = len(self.perm) + 1
        if resized.cols < min_width:
            resized = resized.set_width(min_width)
        results = set()

        # the_perm = self.perm

        # up_perms = [L[0] for L in (ASx(the_perm, self.rows)*ASx(Permutation([]), 1)).keys()]
        # for up_perm in up_perms:

        #     new_bpd = resized.copy()

        # return results
        # Get all crossings in the new row as numpy array for vectorization
        crossings = np.array([(i, j) for (i, j) in resized.all_crossings() if i == self.rows], dtype=np.int32)

        # if len(crossings) == 0:
        #     # No crossings, just rebuild and validate
        #     if resized.is_valid:
        #         results.add(resized.resize(self.rows + 1, column_perm=Permutation([])).set_width(max(resized.rows, len(resized.perm))))
        # else:
        # Generate all 2^n binary masks for subsets
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
            new_bpd = BPD(working_grid, column_perm=resized._column_perm)

            if new_bpd.is_valid and new_bpd.is_reduced:
                results.add(new_bpd.resize(self.rows + 1, column_perm=Permutation([])).set_width(max(new_bpd.rows, len(new_bpd.perm))))

        return results

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

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
from functools import cache
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

    @property
    def is_crossing(self) -> bool:
        """True if this tile is a crossing"""
        return self == TileType.CROSS

    @property
    def is_elbow(self) -> bool:
        """True if this tile is any type of elbow"""
        return self in (TileType.ELBOW_NW, TileType.ELBOW_SE)

    @property
    def is_empty(self) -> bool:
        """True if this tile is empty (pipes go straight)"""
        return self == TileType.BLANK

    @property
    def feeds_right(self) -> bool:
        """True if the horizontal pipe continues to the right"""
        return self in (TileType.HORIZ, TileType.ELBOW_SE, TileType.CROSS)

    @property
    def feeds_up(self) -> bool:
        """True if the vertical pipe continues upwards"""
        return self in (TileType.VERT, TileType.ELBOW_NW, TileType.CROSS)

    @property
    def entrance_from_bottom(self) -> bool:
        """True if a pipe can enter from the bottom"""
        return self in (TileType.VERT, TileType.ELBOW_SE, TileType.CROSS)

    @property
    def entrance_from_left(self) -> bool:
        """True if a pipe can enter from the left"""
        return self in (TileType.HORIZ, TileType.ELBOW_NW, TileType.CROSS)


class BPD(SchubertMonomialGraph, DefaultPrinting):
    """
    Bumpless Pipe Dream representation.

    A bumpless pipe dream is an n×n grid where:
    - TileType.CROSS (1) represents a crossing
    - TileType.BLANK (0) represents an empty box (pipes go straight)
    - For general pipe dreams, can use TileType.ELBOW_* (2-5) for elbows

    Each BPD corresponds to a permutation and has an associated weight.
    """

    def __init__(self, grid, column_perm: Permutation | None = None) -> None:
        """
        Initialize a BPD from a grid.

        Args:
            grid: n×n array-like of TileType values, integers 0-5, or list of lists
        s"""
        self.grid = np.array(grid, dtype=TileType)
        self._column_perm = column_perm if column_perm else Permutation([])
        # Validate grid is square
        # if len(self.grid.shape) != 2 or self.grid.shape[0] != self.grid.shape[1]:
        #     raise ValueError("BPD grid must be square (n×n)")

        self.build()

    @property
    def rows(self) -> int:
        return self.grid.shape[0]

    @property
    def cols(self) -> int:
        return self.grid.shape[1]

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

    @classmethod
    def _get_tbd_tile(cls, left_tile: TileType | None, up_tile: TileType | None) -> TileType:
        if up_tile is None:
            if left_tile is None or not left_tile.feeds_right:
                # if self.DEBUG:
                #     print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
                return TileType.ELBOW_SE
            return TileType.HORIZ
        if left_tile is None:
            if up_tile.entrance_from_bottom:
                return TileType.VERT
            # if self.DEBUG:
            #     print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
            return TileType.ELBOW_SE
        if left_tile.feeds_right:
            if up_tile.entrance_from_bottom:
                # if right_edge:
                #     return TileType.CROSS
                return TileType.ELBOW_NW
            return TileType.HORIZ
        if up_tile.entrance_from_bottom:
            # if right_edge:
            #     return TileType.CROSS
            return TileType.VERT
        # if self.DEBUG:
        #     print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
        return TileType.ELBOW_SE

    DEBUG = False

    def build(self, validate: bool = False) -> None:
        """Build internal structures if needed (currently a placeholder)"""
        for col in range(self.cols):
            for row in range(self.rows):
                if self[row, col] == TileType.TBD:
                    self.grid[row, col] = self._get_tbd_tile(
                        self[row, col - 1] if col > 0 else None,
                        self[row - 1, col] if row > 0 else None,
                    )
                # if self.DEBUG:
                #     print(f"After building cell ({row}, {col}):\n{self}")
        if validate:
            assert self.is_valid(), f"Built BPD is not valid: \n{self}"

    def __len__(self) -> int:
        """Return the size n of the n×n grid"""
        return self.rows

    def __getitem__(self, key) -> TileType | np.ndarray:
        """Access grid elements, casting to TileType"""
        result = self.grid[key]
        # If it's a numpy array (from slicing), cast each element
        if isinstance(result, np.ndarray):
            return result.astype(TileType)
        # Otherwise it's a scalar, cast directly
        return TileType(result)

    def _sympyrepr(self, printer=None) -> str:
        """SymPy repr representation"""
        grid_list = self.grid.tolist()
        if printer is None:
            return f"BPD({grid_list})"
        return f"BPD({printer._print(grid_list)})"

    def _sympystr(self, printer=None) -> str:
        """SymPy str representation of the BPD using tile symbols"""
        # result = []
        # for i in range(self.rows):
        #     row = []
        #     for j in range(self.cols):
        #         tile = TileType(self[i, j])
        #         row.append(str(tile))
        #     result.append("".join(row))
        # if printer is None:
        #     return "\n".join(result)
        return printer._print(pretty(self))

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
        buildperm = []

        for row in range(self.rows):
            # trace from the right
            going_left = True
            current_row = row
            current_col = self.cols - 1
            while current_col >= 0 and current_row < self.rows:
                tile = self[current_row, current_col]
                if tile == TileType.ELBOW_NW:
                    # Elbow NW: go down
                    assert not going_left, "Invalid pipe direction at ELBOW_NW"
                    current_col -= 1
                    going_left = True
                elif tile == TileType.ELBOW_SE:
                    # Elbow SE: go right
                    assert going_left, "Invalid pipe direction at ELBOW_SE"
                    current_row += 1
                    going_left = False
                elif tile == TileType.CROSS:
                    if going_left:
                        current_col -= 1
                    else:
                        current_row += 1
                elif tile == TileType.HORIZ:
                    assert going_left
                    current_col -= 1
                elif tile == TileType.VERT:
                    assert not going_left
                    current_row += 1
                else:
                    raise ValueError(f"Invalid tile type {tile} in BPD")
            buildperm.append(current_col + 1)
        n = max(buildperm, default=1)
        return self._column_perm * Permutation.from_partial(buildperm + [None] * (n - len(buildperm)))

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
        return tuple(int(np.sum(self[i, :] == TileType.BLANK)) for i in range(self.rows))

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

        def is_blank(i, j):
            return corner_sum[i, j] - (corner_sum[i - 1, j - 1] if i > 0 and j > 0 else 0) == 0

        def is_crossing(i, j):
            return corner_sum[i, j] - (corner_sum[i - 1, j - 1] if i > 0 and j > 0 else 0) == 2

        for i in range(n):
            for j in range(n):
                if is_blank(i, j):
                    grid[i, j] = TileType.BLANK
                elif is_crossing(i, j):
                    grid[i, j] = TileType.CROSS
        return cls(grid)

    @classmethod
    def rothe_bpd(cls, perm: Permutation, num_rows: int | None = None) -> BPD:
        if num_rows is None:
            num_rows = len(perm)
        n = max(num_rows, len(perm))
        grid = np.full((num_rows, n), fill_value=TileType.TBD, dtype=TileType)
        bpd = BPD(grid)
        for a, b in perm.diagram:
            bpd.grid[a - 1, b - 1] = TileType.BLANK
        graph = perm.graph
        for i in range(num_rows):
            for j in range(n):
                if bpd[i, j] != TileType.BLANK:
                    if any(tup[0] == i + 1 and tup[1] - 1 < j for tup in graph) and any(tup[0] - 1 < i and tup[1] == j + 1 for tup in graph):
                        bpd.grid[i, j] = TileType.CROSS
        bpd.rebuild()
        return bpd

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
        word = []

        # Read crossings from bottom to top, left to right
        for row in range(self.rows - 1, -1, -1):
            for col in range(self.cols):
                if self[row, col] == 1:
                    # Count pipes weakly northeast of this crossing
                    # Weakly northeast: r <= row and c >= col
                    pipes_northeast = self.cols - col

                    # Word value is pipes_northeast - 1 (1-indexed)
                    word.append(pipes_northeast)

        return tuple(word)

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
            new_grid[:, :width] = self.grid[:, :width]
            bop = BPD(new_grid, column_perm=self._column_perm)
            bop.rebuild()
            assert bop.is_valid, f"Resulting BPD is not valid after reducing width, {pretty(self)} {width} "
            return bop
        if width == self.cols:
            return self
        new_grid = np.full((self.rows, width), fill_value=TileType.TBD, dtype=TileType)
        new_grid[:, : self.cols] = self.grid
        return BPD(new_grid, column_perm=self._column_perm)

    @property
    def is_valid(self) -> bool:
        """
        Check if this is a valid pipe dream.

        A valid pipe dream must satisfy:
        1. Each pipe path from bottom to top does not go out of bounds
        2. The resulting permutation is valid

        Returns:
            True if valid, False otherwise
        """
        # sanity checks
        for row in range(self.rows):
            if self[row, 0] in (TileType.ELBOW_NW, TileType.CROSS, TileType.HORIZ):
                return False
            if self[row, self.cols - 1] in (TileType.ELBOW_NW, TileType.VERT):
                return False
        for col in range(self.cols):
            # if self[self.rows - 1, col] in (TileType.ELBOW_NW, TileType.HORIZ):
            #     return False
            if self[0, col] in (TileType.ELBOW_NW, TileType.VERT, TileType.CROSS):
                return False
        for row in range(1, self.rows - 1):
            for col in range(1, self.cols - 1):
                if self[row, col].feeds_right and not self[row, col + 1].entrance_from_left:
                    return False
                if self[row, col].feeds_up and not self[row - 1, col].entrance_from_bottom:
                    return False
                if self[row, col].entrance_from_left and not self[row, col - 1].feeds_right:
                    return False
                if self[row, col].entrance_from_bottom and not self[row + 1, col].feeds_up:
                    return False
        try:
            self.to_rc_graph()
            return self.perm.inv == sum(self.length_vector)
        except Exception:
            return False

    def __eq__(self, other: object) -> bool:
        """Check equality of two BPDs"""
        if not isinstance(other, BPD):
            return False
        return np.array_equal(self.grid, other.grid)

    def __hash__(self) -> int:
        """Hash for use in sets and dicts"""
        return hash(self.grid.tobytes())

    def copy(self) -> BPD:
        """Create a copy of this BPD"""
        return BPD(self.grid.copy(), column_perm=self._column_perm)

    @property
    def num_crossings(self) -> int:
        """Total number of crossings in the BPD"""
        return int(np.sum(self.grid == TileType.CROSS))

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

    def all_se_elbow_spots(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.ELBOW_SE)

    def all_blank_spots(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.BLANK)

    def all_crossings(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(TileType.CROSS)

    def all_tiles_of_type(self, tile_type: TileType) -> set[tuple[int, int]]:
        tiles = set()
        for i in range(self.rows):
            for j in range(self.cols):
                if self[i, j] == tile_type:
                    tiles.add((i, j))
        return tiles

    def droop_moves(self) -> set[tuple[tuple[int, int], tuple[int, int]]]:
        import itertools

        droop_moves = set()
        for (ri, rj), (bi, bj) in itertools.product(self.all_se_elbow_spots(), self.all_blank_spots()):
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

        D.grid[ri, rj] = TileType.BLANK  # NW-corner
        D.grid[bi, bj] = TileType.ELBOW_NW  # SE-corner
        D.grid[bi, rj] = TileType.ELBOW_SE  # SW-corner
        D.grid[ri, bj] = TileType.ELBOW_SE  # NE-corner

        # top and bottom
        for j in range(rj + 1, bj):
            if self[ri, j] == TileType.HORIZ:
                D.grid[ri, j] = TileType.BLANK
            else:  # self[ri, j] == TileType.CROSS
                D.grid[ri, j] = TileType.VERT
            if self[bi, j] == TileType.BLANK:
                D.grid[bi, j] = TileType.HORIZ
            else:  # self[bi, j] == TileType.VERT
                D.grid[bi, j] = TileType.CROSS
        # left and right
        for i in range(ri + 1, bi):
            if self[i, rj] == TileType.VERT:
                D.grid[i, rj] = TileType.BLANK
            else:  # self[i, rj] == TileType.CROSS
                D.grid[i, rj] = TileType.HORIZ
            if self[i, bj] == TileType.BLANK:
                D.grid[i, bj] = TileType.VERT
            else:  # self[i, bj] == TileType.HORIZ
                D.grid[i, bj] = TileType.CROSS
        D.rebuild()
        return D

    def normalize(self) -> BPD:
        if len(self) >= len(self.perm):
            return self.copy()
        snap_size = max(len(self.perm), self.cols)
        new_grid = np.pad(self.grid, ((0, snap_size - len(self)), (0, max(0, snap_size - self.cols))), constant_values=TileType.TBD)
        bottom_portion = BPD.rothe_bpd(self.perm.min_coset_rep(*(list(range(self.rows)) + list(range(self.rows + 1, snap_size)))), snap_size)
        new_grid[self.rows :, :] = bottom_portion.grid[self.rows :, :]
        ret = BPD(new_grid)
        ret.rebuild()
        return ret

    def pop_op(self) -> tuple[BPD, tuple[int, int]]:
        # --- STEP 0 --- #
        D = self.normalize()
        # check if D has a blank tile (i.e., the coxeter length of D.w is zero)
        if self.perm.inv == 0:
            return D

        # find the first row r with a blank tile
        r = min([i for i in range(D.rows) for j in range(D.cols) if D[i, j] == TileType.BLANK])

        # initialize the mark X at the blank tile (x,y)
        x = r
        y = max([j for i in range(D.rows) for j in range(D.cols) if D[i, j] == TileType.BLANK and i == r])
        while True:
            # --- STEP 1 --- #
            # move the mark to the rightmost blank tile in the block
            j = y
            while j < D.cols and D[x, j] == TileType.BLANK:
                j += 1
            y = j - 1

            # --- STEP 2 --- #
            # find first j-elbow (x_,y+1) with x_>x, if not found then set x_=0 (in which case p==y+1)
            x_ = 0
            for i in range(x + 1, D.rows):
                if D[i, y + 1] == TileType.ELBOW_NW:
                    x_ = i
                    break
            if x_ == 0:  # p==y+1
                break

            for z in range(x + 1, x_):
                if D[z, y] == TileType.BLANK:
                    D.grid[z, y] = TileType.VERT
                    D.grid[z, y + 1] = TileType.BLANK
                elif D[z, y] == TileType.CROSS:
                    continue
                elif D[z, y] == TileType.HORIZ:
                    D.grid[z, y] = TileType.CROSS
                    D.grid[z, y + 1] = TileType.HORIZ
                elif D[z, y] == TileType.VERT:
                    continue
                elif D[z, y] == TileType.ELBOW_SE:
                    D.grid[z, y] = TileType.VERT
                    D.grid[z, y + 1] = TileType.ELBOW_SE
                elif D[z, y] == TileType.ELBOW_NW:
                    D.grid[z, y] = TileType.CROSS
                    D.grid[z, y + 1] = TileType.ELBOW_NW
            D.grid[x, y] = TileType.ELBOW_SE  # NW-corner
            D.grid[x_, y + 1] = TileType.BLANK  # SE-corner
            if D[x_, y] == TileType.ELBOW_SE:  # SW-corner
                D.grid[x_, y] = TileType.VERT
            else:  # D.grid[x_,y] == TileType.HORIZ
                D.grid[x_, y] = TileType.ELBOW_NW
            if D[x, y + 1] == TileType.ELBOW_SE:  # NE-corner
                D.grid[x, y + 1] = TileType.HORIZ
            else:  # D.grid[x,y+1] == TileType.VERT
                D.grid[x, y + 1] = TileType.ELBOW_NW

            # move the mark X to the SE-corner of U
            x = x_
            y = y + 1

        # --- STEP 3 --- #
        a = y  # where (x,y) is the final position of the mark X

        x_ = 0
        for i in range(D.rows - 1, x, -1):
            if D[i, y] == TileType.ELBOW_SE:
                x_ = i
                break

        # copied from above
        for z in range(x + 1, x_):
            if D[z, y] == TileType.BLANK:
                D.grid[z, y] = TileType.VERT
                D.grid[z, y + 1] = TileType.BLANK
            elif D[z, y] == TileType.CROSS:
                continue
            elif D[z, y] == TileType.HORIZ:
                D.grid[z, y] = TileType.CROSS
                D.grid[z, y + 1] = TileType.HORIZ
            elif D[z, y] == TileType.VERT:
                continue
            elif D[z, y] == TileType.ELBOW_SE:
                D.grid[z, y] = TileType.VERT
                D.grid[z, y + 1] = TileType.ELBOW_SE
            elif D[z, y] == TileType.ELBOW_NW:
                D.grid[z, y] = TileType.CROSS
                D.grid[z, y + 1] = TileType.ELBOW_NW

        D.grid[x, y] = TileType.ELBOW_SE  # NW-corner
        D.grid[x_, y + 1] = TileType.ELBOW_SE  # SE-corner
        D.grid[x_, y] = TileType.VERT  # SW-corner
        if D[x, y + 1] == TileType.ELBOW_SE:  # NE-corner
            D.grid[x, y + 1] = TileType.HORIZ
        else:  # D.grid[x,y+1] == TileType.VERT
            D.grid[x, y + 1] = TileType.ELBOW_NW

        D.rebuild()
        # left_simple_ref = Permutation.ref_product(self._column_perm[a])
        # if D.perm != left_simple_ref * self.perm:
        #     D._column_perm = left_simple_ref * self._column_perm
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
            return BPD.from_rc_graph(self.to_rc_graph().resize(new_num_rows), column_perm=column_perm)
        if new_num_rows < len(self.perm):
            return BPD(self.grid[:new_num_rows, :], column_perm=self.column_perm_at_row(new_num_rows - 1) if column_perm is None else column_perm)
        return BPD(self.grid[:new_num_rows, :], column_perm=column_perm)

    @classmethod
    def from_rc_graph(cls, rc_graph, column_perm: Permutation | None = None) -> BPD:
        num_rows = len(rc_graph)
        n = max(num_rows, len(rc_graph.perm))
        bpd = BPD(np.full((n, n), fill_value=TileType.TBD, dtype=TileType))
        coords = [rc_graph.left_to_right_inversion_coords(i) for i in range(rc_graph.perm.inv)]
        coords.reverse()
        # for i, j in coords:
        #     bpd = bpd.inverse_pop_op((i + j - 1, i))
        #     rc_graph = rc_graph.toggle_ref_at(i, j)
        bpd = bpd.inverse_pop_op(*[(i + j - 1, i) for i, j in coords])

        assert bpd.perm.inv == len(bpd.all_blank_spots())
        if len(bpd) != len(rc_graph):
            return bpd.resize(num_rows, column_perm=column_perm)
        return bpd

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

        # Shift the column permutation
        new_column_perm = Permutation([p + shift for p in self._column_perm] if self._column_perm else [])

        return BPD(new_grid, column_perm=new_column_perm)

    def product(self, other: BPD) -> dict[BPD, int]:
        """Compute the product of this BPD with another."""
        from schubmult.utils.perm_utils import add_perm_dict
        other_graph = other.to_rc_graph()
        # other_reduced_compatible = [(a + len(self), r + len(self)) for a, r in other.as_reduced_compatible()]
        # other_reduced_compatible.reverse()
        if self.perm.inv == 0:
            # return {BPD.rothe_bpd(Permutation([]), len(self) + len(other)).inverse_pop_op(*other_reduced_compatible).resize(len(self) + len(other)): 1}
            return {BPD.from_rc_graph(other_graph.prepend(len(self))): 1}
        num_zeros = max(len(other), len(other.perm))
        assert len(self.perm.trimcode) <= len(self), f"{self=}, {self.perm=}"
        base_bpd = self.copy()
        buildup_module = {base_bpd: 1}

        for _ in range(num_zeros):
            new_buildup_module = {}
            for bpd, coeff in buildup_module.items():
                new_buildup_module = add_perm_dict(new_buildup_module, dict.fromkeys(bpd.right_zero_act(), coeff))
            buildup_module = new_buildup_module
        ret_module = {}

        for bpd, coeff in buildup_module.items():
            assert bpd.is_valid, f"Invalid BPD in product buildup: {pretty(bpd)}"
            # try:
            #     new_bpd = bpd.inverse_pop_op(*other_reduced_compatible).resize(len(self) + len(other))
            # except Exception:
            #     continue
            new_rc = RCGraph([*bpd.to_rc_graph()[: len(self)], *other_graph.shiftup(len(self))])
            if new_rc.is_valid:
                new_bpd = BPD.from_rc_graph(new_rc)
                assert len(new_bpd) == len(self) + len(other)

                if new_bpd.is_valid and new_bpd.perm.inv == self.perm.inv + other.perm.inv and len(new_bpd.perm.trimcode) <= len(new_bpd):
                    ret_module = add_perm_dict(ret_module, {new_bpd: coeff})

        return ret_module

    def prod_with_bpd(self, other: BPD) -> BPD:
        """Deprecated: Use product() instead. Returns the single BPD from product dictionary."""
        return self.product(other)

    def inverse_pop_op(self, *interlaced_rc) -> BPD:
        D = self.normalize()

        if len(interlaced_rc) == 2 and isinstance(interlaced_rc[0], int):
            interlaced_rc = [interlaced_rc]
        else:
            interlaced_rc = list(reversed(interlaced_rc))
        while len(interlaced_rc) > 0:
            a, r = interlaced_rc.pop()
            if D.rows <= a or D.rows <= r:
                new_num_rows = max(D.rows, a + 1, r + 1)
                D = D.resize(new_num_rows)
            # find first elbow in column a
            x_ = D.rows - 1
            while x_ >= 0 and D[x_, a] != TileType.ELBOW_SE:
                x_ -= 1
            if x_ == -1:
                raise ValueError("No elbow found in specified column for inverse pop operation when inserting ")
            D.grid[x_, a] = TileType.CROSS
            y = a - 1
            # find x in column y, it will be an SE elbow
            x = x_ - 1
            while x >= 0 and D[x, y] != TileType.ELBOW_SE:
                x -= 1
            if x == -1:
                raise ValueError("No elbow found in specified column for inverse pop operation")
            D.grid[x, y] = TileType.BLANK

            for z in range(x + 1, x_):
                if D[z, y] == TileType.VERT and D[z, y + 1] == TileType.BLANK:
                    D.grid[z, y] = TileType.BLANK
                    D.grid[z, y + 1] = TileType.VERT
                elif D[z, y] == TileType.CROSS and D[z, y + 1] == TileType.ELBOW_NW:
                    D.grid[z, y] = TileType.TBD
                    D.grid[z, y + 1] = TileType.TBD
                elif D[z, y] == TileType.CROSS and D[z, y + 1] == TileType.HORIZ:
                    D.grid[z, y] = TileType.TBD
                    D.grid[z, y + 1] = TileType.CROSS
                elif D[z, y] == TileType.VERT and D[z, y + 1] == TileType.ELBOW_SE:
                    D.grid[z, y] = TileType.ELBOW_SE
                    D.grid[z, y + 1] = TileType.CROSS

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
                while y >= 0 and D[x_, y] == TileType.BLANK:
                    y -= 1
                D.grid[x_, y + 1] = TileType.ELBOW_NW

                # find x at SE elbow
                x = x_ - 1
                while x >= 0 and D[x, y] != TileType.ELBOW_SE:
                    x -= 1
                if x == -1:
                    raise ValueError("No elbow found in specified column for inverse pop operation")

                # [x, y] becomes BLANK
                D.grid[x, y] = TileType.BLANK

                # then this is the fix logic
                for z in range(x + 1, x_):
                    if D[z, y] == TileType.VERT and D[z, y + 1] == TileType.BLANK:
                        D.grid[z, y] = TileType.BLANK
                        D.grid[z, y + 1] = TileType.VERT
                    elif D[z, y] == TileType.CROSS and D[z, y + 1] == TileType.ELBOW_NW:
                        D.grid[z, y] = TileType.TBD
                        D.grid[z, y + 1] = TileType.TBD
                    elif D[z, y] == TileType.CROSS and D[z, y + 1] == TileType.HORIZ:
                        D.grid[z, y] = TileType.TBD
                        D.grid[z, y + 1] = TileType.CROSS
                    elif D[z, y] == TileType.VERT and D[z, y + 1] == TileType.ELBOW_SE:
                        D.grid[z, y] = TileType.ELBOW_SE
                        D.grid[z, y + 1] = TileType.CROSS

                D.rebuild()

        D.rebuild()

        # if D.perm != Permutation.ref_product(self._column_perm[a]) * self.perm:
        #     D._column_perm = Permutation.ref_product(self._column_perm[a]) * self._column_perm

        return D

    def as_reduced_compatible(self):
        work_bpd = self
        ret = []
        while work_bpd.perm.inv > 0:
            work_bpd, (reflection, row) = work_bpd.pop_op()
            ret.append((reflection, row))
        ret.reverse()
        return tuple(ret)

    def rebuild(self) -> None:
        """Rebuild the BPD to resolve any TBD tiles"""
        for i in range(self.rows):
            for j in range(self.cols):
                if self[i, j] not in (TileType.BLANK, TileType.CROSS):
                    self.grid[i, j] = TileType.TBD
        self.build()

    def zero_out_last_row(self) -> BPD:
        return self.resize(self.rows - 1, column_perm=Permutation([]))

    def set_tile(self, i: int, j: int, tile_type: TileType) -> None:
        new_bpd = self.copy()
        new_bpd.grid[i, j] = tile_type
        new_bpd.rebuild()
        return new_bpd


    def right_zero_act(self) -> set[BPD]:
        # find crosses, untransition them
        import itertools
        if not self.is_valid:
            return set()
        resized = self.resize(self.rows + 1)
        if resized.cols < len(self.perm) + 1:
            resized = resized.set_width(len(self.perm) + 1)
        results = set()
        crossings = [(i, j) for (i, j) in resized.all_crossings() if i == self.rows]
        for r in range(len(crossings) + 1):
            for cross_subset in itertools.combinations(crossings, r):
                new_bpd = resized.copy()
                for (i, j) in cross_subset:
                    new_bpd.grid[i, j] = TileType.TBD
                new_bpd.rebuild()
                if new_bpd.is_valid:
                    baggage = new_bpd.resize(self.rows + 1, column_perm=Permutation([])).set_width(max(new_bpd.rows, len(new_bpd.perm)))
                    if baggage.is_valid:
                        results.add(baggage)
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
            return prod(x[i + 1] for i, _ in self.all_blank_spots())
        return prod(x[i + 1] - y[j + 1] for i, j in self.all_blank_spots())

    def to_rc_graph(self) -> RCGraph:
        """
        Convert this BPD to an RC-graph representation.

        Returns:
            RCGraph object (if available in the module)
        """

        # RC-graph rows contain the column positions of crossings in each row
        rows = [[] for _ in range(self.rows)]
        # work_bpd = self
        rc = self.as_reduced_compatible()
        for (reflection, row) in reversed(rc):
            rows[row - 1] = rows[row - 1] + [reflection]

        return RCGraph([tuple(r) for r in rows])

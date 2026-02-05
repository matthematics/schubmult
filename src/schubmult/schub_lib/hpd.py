"""
Bumpless Pipe Dreams (HPD) module
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import IntEnum
from functools import cache, cached_property
from typing import Tuple

import numpy as np
from sympy import pretty
from sympy.printing.defaults import DefaultPrinting

from schubmult.schub_lib.bpd import BPD, TileType
from schubmult.schub_lib.permutation import Permutation
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.schub_lib.schubert_monomial_graph import SchubertMonomialGraph
from schubmult.symbolic import Expr


class HPDTile(IntEnum):
    """
    Enumeration of the possible tile types in a pipe dream.

    Each tile represents how two pipes (horizontal and vertical) interact in a square.
    Note: Whether a tile is "weighty" depends on the row's _id_vector value, not the tile itself.
    """

    TBD = 0  # Placeholder for uninitialized tile
    BLANK = 1  # Both pipes go straight (no crossing, no elbow)
    CROSS = 2  # Pipes cross each other
    HORIZ = 3  # Horizontal pipe
    ELBOW_NW = 4  # Elbow: bottom-right to top-left (╯)
    ELBOW_SE = 5  # Elbow: top-left to bottom-right (╮)
    ELBOW_NE = 6  # Elbow: bottom-right to top-left (╯)
    ELBOW_SW = 7  # Elbow: top-left to bottom-right (╮)
    VERT = 8
    BUMP = 9  # Bump/osculating tile (pipes touch at corner)

    # Legacy aliases for backward compatibility
    UBLANK = 1
    WBLANK = 1
    UCROSS = 2
    WCROSS = 2
    UHORIZ = 3
    WHORIZ = 3

    def __str__(self) -> str:
        """Return base display symbol (use get_display_symbol() for context-aware rendering)"""
        symbols = {
            HPDTile.BLANK: "□",
            HPDTile.CROSS: "┼",
            HPDTile.HORIZ: "─",
            HPDTile.ELBOW_NW: "╯",
            HPDTile.ELBOW_SE: "╭",
            HPDTile.ELBOW_NE: "┕",
            HPDTile.ELBOW_SW: "┑",
            HPDTile.VERT: "│",
            HPDTile.BUMP: "╬",
            HPDTile.TBD: "?",
        }
        return symbols.get(self, "?")

    def get_display_symbol(self, is_weighty: bool) -> str:
        """Get display symbol based on whether the tile is in a weighty row"""
        if is_weighty:
            weighty_symbols = {
                HPDTile.BLANK: "░",
                HPDTile.CROSS: "╋",
                HPDTile.HORIZ: "━",
            }
            return weighty_symbols.get(self, str(self))
        return str(self)

    @classmethod
    def from_tiletype(cls, tile: TileType) -> HPDTile:
        """Convert from TileType to HPDTile"""
        mapping = {
            TileType.BLANK: HPDTile.BLANK,
            TileType.CROSS: HPDTile.CROSS,
            TileType.ELBOW_NW: HPDTile.ELBOW_SW,
            TileType.ELBOW_SE: HPDTile.ELBOW_NE,
            TileType.VERT: HPDTile.VERT,
            TileType.HORIZ: HPDTile.HORIZ,
            TileType.BUMP: HPDTile.BUMP,
        }
        return mapping[tile]

    @cached_property
    def is_crossing(self) -> bool:
        """True if this tile is a crossing"""
        return self == HPDTile.CROSS

    @cached_property
    def is_elbow(self) -> bool:
        """True if this tile is any type of elbow"""
        return self in (HPDTile.ELBOW_NW, HPDTile.ELBOW_SE, HPDTile.ELBOW_NE, HPDTile.ELBOW_SW)

    @cached_property
    def is_empty(self) -> bool:
        """True if this tile is empty (pipes go straight)"""
        return self == HPDTile.BLANK

    # Removed is_weighty - now context-dependent on row's _id_vector
    # Use HPD.is_weighty_position(row, col) instead

    @cached_property
    def feeds_right(self) -> bool:
        """True if the horizontal pipe continues to the right"""
        return self in (HPDTile.HORIZ, HPDTile.ELBOW_SE, HPDTile.CROSS, HPDTile.BUMP, HPDTile.ELBOW_NE)

    @cached_property
    def feeds_up(self) -> bool:
        """True if the vertical pipe continues upwards"""
        return self in (HPDTile.VERT, HPDTile.ELBOW_NW, HPDTile.CROSS, HPDTile.BUMP, HPDTile.ELBOW_NE)

    @cached_property
    def entrance_from_bottom(self) -> bool:
        """True if a pipe can enter from the bottom"""
        return self in (HPDTile.VERT, HPDTile.ELBOW_SE, HPDTile.ELBOW_SW, HPDTile.CROSS, HPDTile.BUMP)

    @cached_property
    def entrance_from_left(self) -> bool:
        """True if a pipe can enter from the left"""
        return self in (HPDTile.HORIZ, HPDTile.ELBOW_NW, HPDTile.CROSS, HPDTile.BUMP, HPDTile.ELBOW_SW)


def _invalidate_grid(grid: np.ndarray) -> None:
    """Helper function to invalidate a grid by setting TBD tiles"""
    # Compute mask once: keep only BLANK and CROSS tiles
    mask = (grid == HPDTile.BLANK) | (grid == HPDTile.CROSS) | (grid == HPDTile.BUMP)
    grid[~mask] = HPDTile.TBD


def _display_grid(grid: np.ndarray) -> str:
    rows = []
    for i in range(grid.shape[0]):
        row_parts = []
        for j in range(grid.shape[1]):
            tile = HPDTile(grid[i, j])
            tile_str = str(tile)
            # Add horizontal extensions for HORIZ tiles
            if tile.entrance_from_left and tile.feeds_right:
                tile_str = "─" + tile_str + "─"
            elif tile.feeds_right:
                tile_str = " " + tile_str + "─"
            elif tile.entrance_from_left:
                tile_str = "─" + tile_str + " "
            else:
                tile_str = " " + tile_str + " "
            row_parts.append(tile_str)
        rows.append("".join(row_parts))
    print("\n".join(rows))


def _classical_bottom_row(weight, length):
    row = np.array([HPDTile.BLANK] * length)
    row[:weight] = HPDTile.HORIZ
    row[weight] = HPDTile.ELBOW_NW
    return row


def _bpd_bottom_row(weight, length):
    row = np.array([HPDTile.BLANK] * length)
    row[weight] = HPDTile.ELBOW_NE
    if weight < length - 1:
        row[weight + 1 :] = HPDTile.HORIZ
    return row


class HPD(SchubertMonomialGraph, DefaultPrinting):
    """
    Bumpless Pipe Dream representation.

    A bumpless pipe dream is an n×n grid where:
    - HPDTile.CROSS (1) represents a crossing
    - HPDTile.BLANK (0) represents an empty box (pipes go straight)
    - For general pipe dreams, can use HPDTile.ELBOW_* (2-5) for elbows

    Each HPD corresponds to a permutation and has an associated weight.
    """

    def __init__(self, grid, id_vector, *, _is_copy=False) -> None:
        """
        Initialize a HPD from a grid.

        Args:
            grid: n×n array-like of HPDTile values, integers 0-5, or list of lists
        s"""
        if _is_copy:
            return
        self._grid = np.array(grid, dtype=HPDTile)
        self._perm = None
        self._valid = None
        self._word = None
        self._unzero_cache = None
        self._id_vector = id_vector
        # self.build()

    def is_classic_row(self, row: int) -> bool:
        """Check if row is a classic row (id_vector[row] == 0)"""
        return self._id_vector[row] == 0

    @classmethod
    def concat(cls, bpd, rc):
        """
        Concatenate a BPD and RCGraph into a HPD.

        Args:
            bpd: BPD instance
            rc: RCGraph instance
        """
        common_size = max(len(bpd.perm), len(rc.perm))
        hpd1 = cls.from_bpd(bpd.resize(common_size))
        hpd2 = cls.from_rc_graph(rc.resize(common_size))
        grid = np.full((bpd.rows + rc.rows, common_size), HPDTile.BLANK, dtype=HPDTile)
        grid[: bpd.rows, :] = hpd1._grid[rc.rows :, :]
        grid[bpd.rows :, :] = hpd2._grid[: rc.rows, :]
        # there is a B C section at bpd.rows - 1, bpd.rows
        br = bpd.rows - 1
        cr = bpd.rows
        for col in range(common_size):
            if grid[br, col] == HPDTile.VERT:
                grid[cr, col] = HPDTile.ELBOW_NW
            elif grid[br, col] == HPDTile.ELBOW_NE:
                grid[cr, col] = HPDTile.HORIZ
            elif grid[br, col] == HPDTile.HORIZ:
                grid[cr, col] = HPDTile.HORIZ
            elif grid[br, col] == HPDTile.BLANK:
                if grid[cr, col] != HPDTile.BLANK:
                    grid[cr, col] = HPDTile.ELBOW_SE
            elif grid[br, col] == HPDTile.ELBOW_SW:
                grid[cr, col] = HPDTile.BUMP
        ret = HPD(grid, (1,) * bpd.rows + (0,) * rc.rows)
        return ret

    def is_weighty_position(self, row: int, col: int) -> bool:
        """Check if a position is in; a weighty row (id_vector[row] == 1)"""
        if self._id_vector[row] == 1:
            return self[row, col] == HPDTile.BLANK
        return self[row, col] == HPDTile.HORIZ or self[row, col] == HPDTile.CROSS

    def row_index_to_label(self, row_index: int) -> int:
        """
        Map physical row index (0-based) to row label.

        Row labels are assigned counterclockwise:
        - _id_vector == 1 rows: labeled 1, 2, ... on RIGHT side, going UPWARD (bottom to top)
        - _id_vector == 0 rows: labeled next, on LEFT side, going DOWNWARD (top to bottom)

        Args:
            row_index: Physical row index (0-based, top to bottom)

        Returns:
            Row label (1-based)
        """
        # Count how many _id_vector == 1 rows exist
        right_rows = [i for i in range(len(self._id_vector)) if self._id_vector[i] == 1]
        left_rows = [i for i in range(len(self._id_vector)) if self._id_vector[i] == 0]

        if self._id_vector[row_index] == 0:
            # Right side: count from bottom to top
            # Find position of row_index in right_rows when reversed

            return left_rows.index(row_index) + 1
            # Left side: continue numbering after right rows, count top to bottom
        offset = len(left_rows)
        position_from_top = len(right_rows) - right_rows.index(row_index)
        return offset + position_from_top

    def row_label_to_index(self, label: int) -> int:
        """
        Map row label (1-based) to physical row index (0-based).

        Inverse of row_index_to_label.

        Args:
            label: Row label (1-based)

        Returns:
            Physical row index (0-based, top to bottom)
        """
        right_rows = [i for i in range(len(self._id_vector)) if self._id_vector[i] == 1]
        left_rows = [i for i in range(len(self._id_vector)) if self._id_vector[i] == 0]

        num_left = len(left_rows)

        if label <= num_left:
            # Right side label: map from bottom to top
            # label 1 is bottommost right row, label num_right is topmost right row
            return left_rows[label - 1]
        # Left side label: map from top to bottom
        # label num_right+1 is topmost left row
        right_position = len(right_rows) - label + num_left
        return right_rows[right_position]

    @classmethod
    def from_bpd(cls, bpd: BPD) -> HPD:
        bpd_grid = np.flipud(bpd.resize(max(bpd.rows, len(bpd.perm)))._grid)
        # Vectorize the conversion function to work on arrays
        vectorized_convert = np.vectorize(lambda tile: HPDTile.from_tiletype(TileType(tile)))
        grid = vectorized_convert(bpd_grid).astype(HPDTile)

        return HPD(grid, (1,) * len(bpd.perm))

    @classmethod
    def from_rc_graph(cls, rc: RCGraph) -> HPD:
        """
        Create a HPD from an RC graph.

        Args:
            rc: RCGraph instance
        """
        n = max(len(rc.perm), rc.rows)
        grid = np.full((n, n), HPDTile.BLANK, dtype=HPDTile)
        for i in range(n):
            for j in range(n - i):
                if rc.has_element(i + 1, j + 1):
                    grid[i, j] = HPDTile.CROSS
                else:
                    grid[i, j] = HPDTile.BUMP
            grid[i, n - 1 - i] = HPDTile.ELBOW_NW
        ret = cls(grid, (0,) * n)
        return ret

    def swap_rows(self, row: int) -> HPD:
        """
        Swap a classic row with the row below it.

        Args:
            row: Index of the classic row to swap (0-based)

        Returns:
            New HPD with the specified rows swapped
        """
        if self._id_vector[row - 1] == self._id_vector[row]:
            raise ValueError("Cannot swap rows, they are the same type")
        if self._id_vector[row - 1] == 1:
            raise ValueError("Top row must be classic (_id_vector == 0) to swap")
        new_grid = self._grid.copy()
        new_id_vector = list(self._id_vector)
        # Swap the rows in the grid
        # new_grid[[row, row + 1], :] = new_grid[[row + 1, row], :]
        # Swap the id_vector entries
        new_id_vector[row - 1], new_id_vector[row] = new_id_vector[row], new_id_vector[row - 1]
        # bpd_row, classic_row = new_grid[row - 1, :], new_grid[row, :]
        # bpd_index, classic_index = row - 1, row
        # if self._id_vector[row - 1] == 0:
        # classic_row, bpd_row = new_grid[row - 1, :], new_grid[row, :]
        classic_index, bpd_index = row - 1, row
        for col in range(self.cols):
            # swap the rows - read from original grid, write to new grid
            classic_tile = self._grid[classic_index, col]
            bpd_tile = self._grid[bpd_index, col]
            if (classic_tile, bpd_tile) == (HPDTile.CROSS, HPDTile.CROSS):
                # nothing to do
                pass
            elif (classic_tile, bpd_tile) == (HPDTile.CROSS, HPDTile.ELBOW_NE):
                new_grid[classic_index, col] = HPDTile.ELBOW_NE
                new_grid[bpd_index, col] = HPDTile.HORIZ
            elif (classic_tile, bpd_tile) == (HPDTile.CROSS, HPDTile.VERT):
                new_grid[classic_index, col] = HPDTile.VERT
                new_grid[bpd_index, col] = HPDTile.CROSS
            elif (classic_tile, bpd_tile) == (HPDTile.HORIZ, HPDTile.HORIZ) and self.pipe_source_labels(classic_index, col)["left"] != self.pipe_source_labels(bpd_index, col)["left"]:
                pass
            elif (classic_tile, bpd_tile) == (HPDTile.HORIZ, HPDTile.HORIZ):
                new_grid[classic_index, col] = HPDTile.BLANK
                new_grid[bpd_index, col] = HPDTile.BLANK
            elif (classic_tile, bpd_tile) == (HPDTile.HORIZ, HPDTile.ELBOW_SW):
                new_grid[classic_index, col] = HPDTile.CROSS
                new_grid[bpd_index, col] = HPDTile.ELBOW_SW
            elif (classic_tile, bpd_tile) == (HPDTile.HORIZ, HPDTile.BLANK):
                new_grid[classic_index, col] = HPDTile.HORIZ
                new_grid[bpd_index, col] = HPDTile.BLANK
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_NW, HPDTile.HORIZ) and self.pipe_source_labels(classic_index, col)["left"] != self.pipe_source_labels(bpd_index, col)[
                "left"
            ]:  # trace pipes
                new_grid[classic_index, col] = HPDTile.CROSS
                new_grid[bpd_index, col] = HPDTile.ELBOW_NW
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_NW, HPDTile.HORIZ):
                new_grid[classic_index, col] = HPDTile.ELBOW_NE
                new_grid[bpd_index, col] = HPDTile.BLANK
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_NW, HPDTile.ELBOW_SW) and self.pipe_source_labels(classic_index, col)["left"] != self.pipe_source_labels(bpd_index, col)[
                "left"
            ]:  # trace pipes
                new_grid[classic_index, col] = HPDTile.CROSS
                new_grid[bpd_index, col] = HPDTile.BUMP
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_NW, HPDTile.ELBOW_SW):
                new_grid[classic_index, col] = HPDTile.ELBOW_NE
                new_grid[bpd_index, col] = HPDTile.ELBOW_SE
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_NW, HPDTile.BLANK):
                new_grid[classic_index, col] = HPDTile.ELBOW_NE
                new_grid[bpd_index, col] = HPDTile.HORIZ
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_SE, HPDTile.CROSS):
                new_grid[classic_index, col] = HPDTile.HORIZ
                new_grid[bpd_index, col] = HPDTile.ELBOW_SE
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_SE, HPDTile.ELBOW_NE):
                new_grid[classic_index, col] = HPDTile.ELBOW_SW
                new_grid[bpd_index, col] = HPDTile.ELBOW_NW
            elif (classic_tile, bpd_tile) == (HPDTile.ELBOW_SE, HPDTile.VERT):
                new_grid[classic_index, col] = HPDTile.ELBOW_SW
                new_grid[bpd_index, col] = HPDTile.BUMP
            elif (
                (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.CROSS)
                and self.pipe_source_labels(classic_index, col)["left"] != self.pipe_source_labels(classic_index, col)["right"]
                and self.pipe_source_labels(bpd_index, col)["top"] != self.pipe_source_labels(classic_index, col)["bottom"] != self.pipe_source_labels(bpd_index, col)["left"]
            ):  # trace pipes
                new_grid[classic_index, col] = HPDTile.CROSS
                new_grid[bpd_index, col] = HPDTile.BUMP
            elif (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.CROSS):
                new_grid[classic_index, col] = HPDTile.ELBOW_NE
                new_grid[bpd_index, col] = HPDTile.ELBOW_SE
            elif (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.ELBOW_NE):
                new_grid[classic_index, col] = HPDTile.VERT
                new_grid[bpd_index, col] = HPDTile.ELBOW_NW
            elif (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.VERT):
                new_grid[classic_index, col] = HPDTile.VERT
                new_grid[bpd_index, col] = HPDTile.BUMP
            elif (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.ELBOW_SE):
                new_grid[classic_index, col] = HPDTile.ELBOW_SW
                new_grid[bpd_index, col] = HPDTile.VERT
            elif (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.HORIZ):
                new_grid[classic_index, col] = HPDTile.HORIZ
                new_grid[bpd_index, col] = HPDTile.VERT
            elif (classic_tile, bpd_tile) == (HPDTile.BUMP, HPDTile.BLANK):
                new_grid[classic_index, col] = HPDTile.BLANK
                new_grid[bpd_index, col] = HPDTile.HORIZ
            elif (classic_tile, bpd_tile) == (HPDTile.BLANK, HPDTile.HORIZ):
                new_grid[classic_index, col] = HPDTile.HORIZ
                new_grid[bpd_index, col] = HPDTile.BLANK
            elif (classic_tile, bpd_tile) == (HPDTile.BLANK, HPDTile.ELBOW_SW):
                new_grid[classic_index, col] = HPDTile.HORIZ
                new_grid[bpd_index, col] = HPDTile.ELBOW_SE
            elif (classic_tile, bpd_tile) == (HPDTile.BLANK, HPDTile.BLANK):
                new_grid[classic_index, col] = HPDTile.HORIZ
                new_grid[bpd_index, col] = HPDTile.HORIZ
            else:
                raise ValueError(f"Cannot swap rows at column {col}: ({classic_tile}, {bpd_tile})")
            # ROW 4
        # _display_grid(new_grid)
        ret = HPD(new_grid, tuple(new_id_vector))
        return ret

    def pipe_source_labels(self, row: int, col: int) -> dict[str, int | None]:
        """
        Determine which row label(s) the pipe(s) at position (row, col) came from.

        Returns a dict with keys 'top', 'bottom', 'left', 'right' indicating which
        row label the pipe on each side of the tile belongs to (None if no pipe on that side).

        Invariants:
        - For non-BUMP/CROSS tiles: exactly 2 non-None sides with equal values
        - For BUMP: right == bottom, left == top
        - For CROSS: top == bottom, left == right

        Args:
            row: Physical row index (0-based)
            col: Physical column index (0-based)

        Returns:
            Dict with keys 'top', 'bottom', 'left', 'right' mapping to row labels or None
        """

        def trace_pipe_to_entry(start_row: int, start_col: int, entry_side: str) -> int | None:
            """
            Trace a pipe from a given tile side to its entry point.

            entry_side: which side of the starting tile the pipe enters from ('left', 'right', 'top', 'bottom')
            Returns: row label where pipe enters, or None if can't trace
            """
            # Start by moving into the tile from the entry side
            curr_row, curr_col = start_row, start_col
            came_from = entry_side  # Which side we entered the current tile from

            visited = set()

            while True:
                # Check for infinite loop
                if (curr_row, curr_col, came_from) in visited:
                    return None
                visited.add((curr_row, curr_col, came_from))

                # Check if out of bounds
                if not (0 <= curr_row < self.rows and 0 <= curr_col < self.cols):
                    return None

                tile = HPDTile(self[curr_row, curr_col])

                # Determine exit side based on tile type and entry side
                if tile == HPDTile.BLANK:
                    return None  # No connection

                if tile == HPDTile.HORIZ:
                    if came_from == "right":
                        exit_side = "left"
                    elif came_from == "left":
                        exit_side = "right"
                    else:
                        return None  # Inconsistent

                elif tile == HPDTile.VERT:
                    if came_from == "top":
                        exit_side = "bottom"
                    elif came_from == "bottom":
                        exit_side = "top"
                    else:
                        return None  # Inconsistent

                elif tile == HPDTile.CROSS:
                    if came_from == "left":
                        exit_side = "right"
                    elif came_from == "right":
                        exit_side = "left"
                    elif came_from == "top":
                        exit_side = "bottom"
                    elif came_from == "bottom":
                        exit_side = "top"
                    else:
                        return None

                elif tile == HPDTile.BUMP:
                    if came_from == "left":
                        exit_side = "top"
                    elif came_from == "top":
                        exit_side = "left"
                    elif came_from == "right":
                        exit_side = "bottom"
                    elif came_from == "bottom":
                        exit_side = "right"
                    else:
                        return None

                elif tile == HPDTile.ELBOW_NW:
                    # Connects left and top
                    if came_from == "left":
                        exit_side = "top"
                    elif came_from == "top":
                        exit_side = "left"
                    else:
                        return None

                elif tile == HPDTile.ELBOW_SE:
                    # Connects bottom and right
                    if came_from == "bottom":
                        exit_side = "right"
                    elif came_from == "right":
                        exit_side = "bottom"
                    else:
                        return None

                elif tile == HPDTile.ELBOW_NE:
                    # Connects right and top
                    if came_from == "right":
                        exit_side = "top"
                    elif came_from == "top":
                        exit_side = "right"
                    else:
                        return None

                elif tile == HPDTile.ELBOW_SW:
                    # Connects bottom and left
                    if came_from == "bottom":
                        exit_side = "left"
                    elif came_from == "left":
                        exit_side = "bottom"
                    else:
                        return None

                else:
                    return None  # Unknown tile

                # Move to next tile based on exit_side
                next_row, next_col = curr_row, curr_col
                # Move to next tile based on exit_side
                next_row, next_col = curr_row, curr_col

                if exit_side == "left":
                    next_col -= 1
                elif exit_side == "right":
                    next_col += 1
                elif exit_side == "top":
                    next_row -= 1
                elif exit_side == "bottom":
                    next_row += 1

                # Check if we've exited the grid (found entry point)
                if next_col < 0:
                    # Exited left edge
                    if self._id_vector[curr_row] == 0:  # Classic row entry
                        return self.row_index_to_label(curr_row)
                    return None  # Not an entry point
                if next_col >= self.cols:
                    # Exited right edge
                    if self._id_vector[curr_row] == 1:  # BPD row entry
                        return self.row_index_to_label(curr_row)
                    return None  # Not an entry point
                if next_row < 0:
                    # Exited top edge - shouldn't happen
                    return None
                if next_row >= self.rows:
                    # Exited bottom edge - bottom entry point
                    return self.row_index_to_label(curr_row)

                # Continue to next tile
                curr_row, curr_col = next_row, next_col
                if exit_side == "left":
                    came_from = "right"
                elif exit_side == "right":
                    came_from = "left"
                elif exit_side == "top":
                    came_from = "bottom"
                elif exit_side == "bottom":
                    came_from = "top"

        tile = HPDTile(self[row, col])
        result = {"top": None, "bottom": None, "left": None, "right": None}

        # Trace each side based on tile connectivity
        if tile == HPDTile.BLANK:
            pass

        elif tile == HPDTile.HORIZ:
            label = trace_pipe_to_entry(row, col, "left")
            if label is None:
                label = trace_pipe_to_entry(row, col, "right")
            result["left"] = label
            result["right"] = label

        elif tile == HPDTile.VERT:
            label = trace_pipe_to_entry(row, col, "bottom")
            if label is None:
                label = trace_pipe_to_entry(row, col, "top")
            result["bottom"] = label
            result["top"] = label

        elif tile == HPDTile.ELBOW_NW:
            # Connects left and top
            label = trace_pipe_to_entry(row, col, "left")
            if label is None:
                label = trace_pipe_to_entry(row, col, "top")
            result["left"] = label
            result["top"] = label

        elif tile == HPDTile.ELBOW_SE:
            # Connects bottom and right
            label = trace_pipe_to_entry(row, col, "bottom")
            if label is None:
                label = trace_pipe_to_entry(row, col, "right")
            result["bottom"] = label
            result["right"] = label

        elif tile == HPDTile.ELBOW_NE:
            # Connects right and top
            label = trace_pipe_to_entry(row, col, "right")
            if label is None:
                label = trace_pipe_to_entry(row, col, "top")
            result["right"] = label
            result["top"] = label

        elif tile == HPDTile.ELBOW_SW:
            # Connects bottom and left
            label = trace_pipe_to_entry(row, col, "bottom")
            if label is None:
                label = trace_pipe_to_entry(row, col, "left")
            result["bottom"] = label
            result["left"] = label

        elif tile == HPDTile.CROSS:
            vert_label = trace_pipe_to_entry(row, col, "bottom")
            if vert_label is None:
                vert_label = trace_pipe_to_entry(row, col, "top")
            horiz_label = trace_pipe_to_entry(row, col, "left")
            if horiz_label is None:
                horiz_label = trace_pipe_to_entry(row, col, "right")
            result["bottom"] = vert_label
            result["top"] = vert_label
            result["left"] = horiz_label
            result["right"] = horiz_label

        elif tile == HPDTile.BUMP:
            right_bottom_label = trace_pipe_to_entry(row, col, "right")
            if right_bottom_label is None:
                right_bottom_label = trace_pipe_to_entry(row, col, "bottom")
            left_top_label = trace_pipe_to_entry(row, col, "left")
            if left_top_label is None:
                left_top_label = trace_pipe_to_entry(row, col, "top")
            result["right"] = right_bottom_label
            result["bottom"] = right_bottom_label
            result["left"] = left_top_label
            result["top"] = left_top_label

        return result

    def _invalidate_cache(self):
        self._perm = None
        self._valid = None
        self._word = None
        self._unzero_cache = None

    @property
    def rows(self) -> int:
        return self._grid.shape[0]

    @property
    def cols(self) -> int:
        return self._grid.shape[1]

    @classmethod
    @cache
    def all_bpds(cls, w: Permutation, length: int | None = None) -> set[HPD]:
        if length is None:
            length = len(w)
        pipes = set()
        new_pipes = [HPD.rothe_bpd(w, length)]

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

    def delete_row(self, row: int) -> HPD:
        new_grid = self._grid.copy()
        # bad_crosses = sorted([(i, j) for (i, j) in self.all_crossings() if 1 in set(self.right_root_at(i, j))])
        # for (i, j) in bad_crosses:
        #     new_grid[i, j: self.cols - 1] = new_grid[i, j + 1 :]
        col = self.cols - 1
        current_row = row - 1

        # if current_row > 0:
        #     for check_col in range(self.cols - 1):
        #         if self[current_row, check_col] == HPDTile.ELBOW_NW:
        #             if self[current_row - 1, check_col] == HPDTile.CROSS:
        #                 new_grid[current_row - 1, check_col] = HPDTile.HORIZ
        #             elif self[current_row - 1, check_col] == HPDTile.VERT:
        #                 new_grid[current_row - 1, check_col] = HPDTile.ELBOW_NW
        #             elif self[current_row - 1, check_col] == HPDTile.ELBOW_SE:
        #                 new_grid[current_row - 1, check_col] = HPDTile.HORIZ
        #             #new_grid[current_row - 1, check_col] = HPDTile.VERT

        going_left = True
        # if current_row < self.rows - 1:
        #     for check_col in range(self.cols - 1):
        #         if self[current_row + 1, check_col] == HPDTile.CROSS:
        #             new_grid[current_row + 1, check_col] = HPDTile.HORIZ
        #         elif self[current_row + 1, check_col] == HPDTile.VERT:
        #             new_grid[current_row + 1, check_col] = HPDTile.ELBOW_SE
        #         elif self[current_row + 1, check_col] == HPDTile.ELBOW_NW:
        #             new_grid[current_row + 1, check_col] = HPDTile.HORIZ
        # if current_row > 0:
        #     for check_col in range(self.cols - 1):
        #         if self[current_row - 1, check_col] == HPDTile.CROSS:
        #             new_grid[current_row - 1, check_col] = HPDTile.HORIZ
        #         elif self[current_row - 1, check_col] == HPDTile.VERT:
        #             new_grid[current_row - 1, check_col] = HPDTile.ELBOW_NW
        #         elif self[current_row - 1, check_col] == HPDTile.ELBOW_SE:
        #             new_grid[current_row - 1, check_col] = HPDTile.HORIZ
        # while current_row > 0 and new_grid[current_row, col] == HPDTile.CROSS:
        #     current_row -= 1
        while current_row < self.rows:
            #     new_grid[current_row, col : self.cols - 1] = new_grid[current_row, col + 1 :]
            # if deleted_tile.entrance_from_left:
            # for j in range(self.cols):
            #     if j > change_col - 1:
            #         grid[n - 1 - row, j] = HPDTile.CROSS
            # for check_row in range(self.rows):
            #     if check_row == row - 1:
            #         continue
            #     if col > 0:
            #         if self[check_row, col - 1] == HPDTile.HORIZ:
            #             new_grid[check_row, col - 1] = HPDTile.ELBOW_NW
            #         elif self[check_row, col - 1] == HPDTile.ELBOW_SE:
            #             new_grid[check_row, col - 1] = HPDTile.VERT
            #         elif self[check_row, col - 1] == HPDTile.CROSS:
            #             new_grid[check_row, col - 1] = HPDTile.VERT
            #     # else:
            #     #     if self[check_row, col] == HPDTile.HORIZ:
            #     #         new_grid[check_row, col] = HPDTile.ELBOW_SE
            #     #     elif self[check_row, col] == HPDTile.ELBOW_NW:
            #     #         new_grid[check_row, col] = HPDTile.VERT
            #     #     elif self[check_row, col] == HPDTile.CROSS:
            #     #         new_grid[check_row, col] = HPDTile.VERT
            #     if col < self.cols - 1:
            #         if self[check_row, col] == HPDTile.HORIZ:
            #             new_grid[check_row, col] = HPDTile.ELBOW_SE
            #         elif self[check_row, col] == HPDTile.ELBOW_NW:
            #             new_grid[check_row, col] = HPDTile.VERT
            #         elif self[check_row, col] == HPDTile.CROSS:
            #             new_grid[check_row, col] = HPDTile.VERT
            #     else:
            #         if self[check_row, col] == HPDTile.HORIZ:
            #             new_grid[check_row, col] = HPDTile.ELBOW_SE
            #         elif self[check_row, col] == HPDTile.CROSS:
            #             new_grid[check_row, col] = HPDTile.HORIZ
            new_grid[current_row, col] = HPDTile.TBD
            if self[current_row, col] == HPDTile.HORIZ:
                col -= 1
                going_left = True
            elif self[current_row, col] == HPDTile.VERT:
                current_row += 1
                going_left = False
            elif self[current_row, col] == HPDTile.CROSS:
                if going_left:
                    col -= 1

                else:
                    current_row += 1
            elif self[current_row, col] == HPDTile.ELBOW_SE:
                current_row += 1
                going_left = False
            elif self[current_row, col] == HPDTile.ELBOW_NW:
                col -= 1
                going_left = True
        # new_grid = np.concatenate([new_grid[:row - 1, : self.cols - 1], new_grid[row:, : self.cols - 1]], axis=0)
        new_grid = np.concatenate([new_grid[: row - 1, : self.cols - 1], new_grid[row:, : self.cols - 1]], axis=0)
        # new_grid = np.delete(new_grid, row - 1, axis=0)

        ret = HPD(new_grid)
        ret.rebuild()
        return ret

    # def interlace(self, other: HPD, start_row: int) -> HPD:
    #     cut_list = [self.resize(i).perm for i in range(start_row)]

    #     for i in range(other.rows - 1):
    #         other_cut0 = other.resize(i)
    #         other_cut = other.resize(i + 1)
    #         #cut_list.append(cut_list[-1] * (~other_cut0.perm * other_cut.perm))
    #         cut_list
    #     the_perm = cut_list[-1]
    #     n = len(the_perm)
    #     while len(cut_list) < n:
    #         cut_list.append(cut_list[-1])
    #     pathfat = tuple([cut_list[i] * Permutation.w0(n - i).shiftup(i) for i in range(n - 1, -1, -1)])
    #     # print(pathfat)
    #     return HPD.from_bruhat_path(pathfat).resize(start_row + other.rows)
    # perm = self.perm * other.perm.shiftup(start_row)
    # solf = self.resize(len(perm))
    # new_grid = solf._grid.copy()

    # last_tbd_col = 0
    # first_tbd_col = -1
    # for col in range(new_grid.shape[1]):
    #     if HPDTile(new_grid[start_row - 1, col]).entrance_from_bottom:
    #         new_grid[start_row:, col] = HPDTile.VERT
    #     else:
    #         new_grid[start_row:, col] = HPDTile.TBD
    #         if first_tbd_col == -1:
    #             first_tbd_col = col
    #         last_tbd_col = col
    # _display_grid(new_grid)
    # for row in range(other.rows):
    #     other_col = other.cols - 1
    #     for col in range(new_grid.shape[1] - 1, first_tbd_col - 1, -1):
    #         if col > last_tbd_col:
    #             new_grid[row + start_row, col] = HPDTile.CROSS
    #             continue

    # if other_col >= other.cols:
    #     if HPDTile(new_grid[row, col]) == HPDTile.TBD:
    #         new_grid[row + self.rows, col] = HPDTile.HORIZ
    #     elif HPDTile(new_grid[row, col]) == HPDTile.VERT:
    #         new_grid[row + self.rows, col] = HPDTile.CROSS
    #     continue
    # wasafeed = False
    # if other_col == 0:
    #     while col < new_grid.shape[1] and HPDTile(new_grid[row + self.rows, col]) == HPDTile.TBD:
    #         new_grid[row + self.rows, col] = HPDTile.HORIZ
    #         col += 1
    #         wasafeed = True
    # while col < new_grid.shape[1] and HPDTile(new_grid[row + self.rows, col]) != HPDTile.TBD:
    #     if (wasafeed and other_col == 0) or (other_col > 0 and other[row, other_col].entrance_from_left  and other[row, other_col - 1].feeds_right):
    #         new_grid[row + start_row, col] = HPDTile.CROSS
    #     col += 1
    #         if HPDTile(new_grid[row + start_row, col]) == HPDTile.TBD:
    #             new_grid[row + start_row, col] = other[row, other_col]
    #             other_col -= 1
    #         elif other[row, other_col].feeds_right:
    #             new_grid[row + start_row, col] = HPDTile.CROSS
    #         _display_grid(new_grid)
    # return HPD(new_grid)

    def append(self, other: HPD) -> HPD:
        perm = self.perm * other.perm.shiftup(self.rows)
        new_grid = np.full((self.rows + other.rows, max(self.cols, len(perm))), HPDTile.TBD, dtype=HPDTile)
        solf = self.resize(new_grid.shape[1])
        new_grid[: self.rows, :] = solf._grid[: self.rows, :]
        # new_grid[self.rows :, self.cols :] = other._grid
        last_tbd_col = 0
        first_tbd_col = -1
        for col in range(new_grid.shape[1]):
            if HPDTile(new_grid[self.rows - 1, col]).entrance_from_bottom:
                new_grid[self.rows :, col] = HPDTile.VERT
            else:
                if first_tbd_col == -1:
                    first_tbd_col = col
                last_tbd_col = col
        _display_grid(new_grid)
        for row in range(other.rows):
            other_col = other.cols - 1
            for col in range(new_grid.shape[1] - 1, first_tbd_col - 1, -1):
                if col > last_tbd_col:
                    new_grid[row + self.rows, col] = HPDTile.CROSS
                    continue

                # if other_col >= other.cols:
                #     if HPDTile(new_grid[row, col]) == HPDTile.TBD:
                #         new_grid[row + self.rows, col] = HPDTile.HORIZ
                #     elif HPDTile(new_grid[row, col]) == HPDTile.VERT:
                #         new_grid[row + self.rows, col] = HPDTile.CROSS
                #     continue
                # wasafeed = False
                # if other_col == 0:
                #     while col < new_grid.shape[1] and HPDTile(new_grid[row + self.rows, col]) == HPDTile.TBD:
                #         new_grid[row + self.rows, col] = HPDTile.HORIZ
                #         col += 1
                #         wasafeed = True
                # while col < new_grid.shape[1] and HPDTile(new_grid[row + self.rows, col]) != HPDTile.TBD:
                #     if (wasafeed and other_col == 0) or (other_col > 0 and other[row, other_col].entrance_from_left  and other[row, other_col - 1].feeds_right):
                #         new_grid[row + self.rows, col] = HPDTile.CROSS
                #     col += 1
                if HPDTile(new_grid[row + self.rows, col]) == HPDTile.TBD:
                    new_grid[row + self.rows, col] = other[row, other_col]
                    other_col -= 1
                elif other[row, other_col].feeds_right:
                    new_grid[row + self.rows, col] = HPDTile.CROSS
                _display_grid(new_grid)
        return HPD(new_grid)

    @classmethod
    def from_bruhat_path(cls, path: Sequence[Permutation]) -> HPD:
        """
        Create a HPD from a Bruhat path.
        """
        n = len(path)
        grid = np.full((n, n), HPDTile.TBD, dtype=HPDTile)

        for row in range(len(path) - 1, 0, -1):
            grid[n - 1 - row, :] = HPD.row_from_k_chain(path[row - 1], path[row], n - row, n)
            change_col = path[row][n - row - 1]
            if grid[n - 1 - row, change_col - 1] == HPDTile.ELBOW_NW:
                grid[n - 1 - row, change_col - 1] = HPDTile.HORIZ
            elif grid[n - 1 - row, change_col - 1] == HPDTile.VERT:
                grid[n - 1 - row, change_col - 1] = HPDTile.ELBOW_SE
            for j in range(n):
                if j > change_col - 1:
                    grid[n - 1 - row, j] = HPDTile.CROSS
        change_col = path[0][n - 1]
        for j in range(n):
            if j > change_col - 1:
                grid[n - 1, j] = HPDTile.CROSS
            # change (1, w(k)) from
            # col = permo[rcol] - 1
            # if col <= rowand path[row][col] == path[row - 1][col]:
            #     grid[row - 1, col] = HPDTile.BLANK
            # elif col > row and path[row][col] == path[row - 1][col] and path[row][col] != Permutation.w0(n)[col]:
            #     grid[row - 1, col] = HPDTile.CROSS
        ret = cls(grid)
        # ret.rebuild()
        return ret

    @staticmethod
    def row_from_k_chain(u: Permutation, w: Permutation, k: int, n: int) -> np.ndarray:
        """
        Construct a single row of tiles from a k-chain according to Definition 3.15.

        Given two permutations u and w where u ≤ w in Bruhat order, finds a maximal
        k-chain from u to w and constructs a row of n tiles based on that chain.

        Args:
            u: Starting permutation
            w: Target permutation (must satisfy u ≤ w in Bruhat order)
            k: The chain parameter (k >= 1)

        Returns:
            1D numpy array of HPDTile values representing the row

        Definition 3.15 cases (for tile at position (row, c)):
        - If chain swaps c with larger but not smaller: ELBOW_SE (⌜)
        - If chain swaps c with both larger and smaller: CROSS (╋)
        - If chain swaps c with smaller but not larger: ELBOW_NW (⌟)
        - If c not among first k numbers of w: BLANK (□)
        - If chain swaps values a,b with a < c < b: BUMP (╬)
        - Otherwise: CROSS (■)
        """
        row_tiles = np.full(n, HPDTile.TBD, dtype=HPDTile)

        # Find a maximal k-chain from u to w
        # A k-chain only swaps among the first k positions

        first_k_numbers = {w[i] for i in range(k)}

        cycles = ((~u) * w).get_cycles()

        bacon = u
        swaps_by_a = {}
        swaps_by_b = {}
        for cyc in cycles:
            cyc = [*cyc]
            a = min(cyc)
            cyc.remove(a)
            while len(cyc) > 0:
                bindex = min(list(range(len(cyc))), key=lambda x: bacon[cyc[x] - 1])
                b = cyc.pop(bindex)
                assert bacon[a - 1] < bacon[b - 1]
                swaps_by_a[bacon[a - 1]] = bacon[b - 1]
                swaps_by_b[bacon[b - 1]] = bacon[a - 1]
                bacon = bacon.swap(a - 1, b - 1)
        assert bacon == w
        # For each column c (1-indexed in definition, 0-indexed in array)
        for c in range(1, n + 1):
            # Determine what values the chain swaps with c

            swaps_with_larger = c in swaps_by_a
            swaps_with_smaller = c in swaps_by_b

            # Determine tile type based on Definition 3.15
            tile = HPDTile.TBD

            # Chain swaps c
            if swaps_with_larger and not swaps_with_smaller:
                tile = HPDTile.ELBOW_SE if c in first_k_numbers else HPDTile.ELBOW_NE
            elif swaps_with_smaller and not swaps_with_larger:
                tile = HPDTile.ELBOW_NW if c in first_k_numbers else HPDTile.ELBOW_SW
            elif swaps_with_larger and swaps_with_smaller:
                tile = HPDTile.HORIZ
            else:
                # c stays unswapped
                if c not in first_k_numbers:
                    tile = HPDTile.BLANK
                else:
                    # Check if chain ever swaps values a,b with a < c < b
                    swaps_bracketing = any(swaps_by_a[a] > c > a for a in swaps_by_a)

                    if swaps_bracketing:
                        tile = HPDTile.CROSS
                    else:
                        tile = HPDTile.BLANK

            # Set the tile at column c-1 (0-indexed)
            row_tiles[c - 1] = tile

        return row_tiles

    # @staticmethod
    # def _find_k_chain(u: Permutation, w: Permutation, k: int) -> list[Permutation]:
    #     """
    #     Find a maximal k-chain from u to w.

    #     A k-chain is a saturated chain in Bruhat order where all transpositions
    #     involve at least one position among the first k positions.

    #     Args:
    #         u: Starting permutation
    #         w: Ending permutation (must satisfy u ≤ w in Bruhat order)
    #         k: Chain parameter (only swap positions involving first k positions)

    #     Returns:
    #         List of permutations forming the k-chain from u to w
    #     """
    #     if u == w:
    #         return [u]

    #     chain = [u]
    #     current = u

    #     # Greedily build a chain by finding valid covers
    #     while current != w:
    #         found_next = False
    #         # Try all possible transpositions that could move us closer to w
    #         for i in range(len(current)):
    #             for j in range(i + 1, len(current) + 1):
    #                 # Check if this is a valid k-chain move (involves first k positions)
    #                 if i >= k and j >= k:
    #                     continue

    #                 # Try swapping positions i and j
    #                 candidate = current.swap(i, j)

    #                 # Check if this is a Bruhat cover and moves us toward w
    #                 if candidate.bruhat_leq(w) and current.inv + 1 == candidate.inv:
    #                     current = candidate
    #                     chain.append(current)
    #                     found_next = True
    #                     break
    #             if found_next:
    #                 break

    #         if not found_next:
    #             # No valid move found; this shouldn't happen if u ≤ w
    #             raise ValueError(f"Cannot find k-chain from {u} to {w} with k={k}")

    #     return chain

    # @classmethod
    # def from_k_chains(cls, w: Permutation, u_dict: dict[int, Permutation]) -> HPD:
    #     """
    #     Create a HPD from a collection of starting permutations, one for each row.

    #     For row k (1-indexed), construct a k-chain from u_k to w and use it
    #     to determine the tiles in that row.

    #     Args:
    #         w: The target permutation
    #         u_dict: Dictionary mapping k (row number, 1-indexed) to starting permutation u_k

    #     Returns:
    #         HPD object constructed from the k-chains
    #     """
    #     n = len(w)
    #     grid = np.full((n, n), HPDTile.TBD, dtype=HPDTile)

    #     for k in range(1, n + 1):
    #         if k not in u_dict:
    #             raise ValueError(f"Missing starting permutation for row {k}")
    #         u_k = u_dict[k]
    #         # Construct row k-1 (0-indexed) using a k-chain from u_k to w
    #         grid[k - 1, :] = cls.from_k_chain(u_k, w, k)

    #     ret = cls(grid)
    #     ret.rebuild()
    #     return ret

    def to_bruhat_path(self):
        n = len(self.perm)
        bigself = self.resize(n)
        return tuple([bigself.resize(i).perm * Permutation.w0(n - i).shiftup(i) for i in range(n - 1, -1, -1)])

    # def cheat_delete_row(self, row: int) -> HPD:
    #     resize_array = []
    #     for i in range(row):
    #         resize_array.append(self.resize(i))
    #     rc = self.to_rc_graph().rowrange(row).shiftup(row - 1)
    #     return ret

    @classmethod
    def _init_tbd_lookup(cls):
        """Initialize the static TBD tile lookup table."""
        if cls._TBD_LOOKUP is not None:
            return

        # Create lookup table: shape (9, 9) for tile types 0-7 (including BUMP=7) plus None (-1 -> index 8)
        lookup = np.full((9, 9), HPDTile.ELBOW_SE, dtype=HPDTile)

        # Map None to index 8
        none_idx = 8

        # For each combination of (left_tile, up_tile), compute result
        # Range includes -1 (None) through 7 (BUMP)
        for left_val in range(-1, 8):  # -1 for None, 0-7 for tile types (TBD, BLANK, CROSS, HORIZ, ELBOW_NW, ELBOW_SE, VERT, BUMP)
            left_idx = none_idx if left_val == -1 else left_val
            left_feeds_right = False if left_val == -1 else HPDTile(left_val).feeds_right

            for up_val in range(-1, 8):
                up_idx = none_idx if up_val == -1 else up_val
                up_entrance = False if up_val == -1 else HPDTile(up_val).entrance_from_bottom

                # Apply the logic from _get_tbd_tile
                if up_val == -1:  # up_tile is None
                    if left_val == -1 or not left_feeds_right:
                        lookup[left_idx, up_idx] = HPDTile.ELBOW_SE
                    else:
                        lookup[left_idx, up_idx] = HPDTile.HORIZ
                elif left_val == -1:  # left_tile is None, up_tile is not
                    if up_entrance:
                        lookup[left_idx, up_idx] = HPDTile.VERT
                    else:
                        lookup[left_idx, up_idx] = HPDTile.ELBOW_SE
                else:  # Both not None
                    if left_feeds_right:
                        if up_entrance:
                            lookup[left_idx, up_idx] = HPDTile.ELBOW_NW
                        else:
                            lookup[left_idx, up_idx] = HPDTile.HORIZ
                    else:
                        if up_entrance:
                            lookup[left_idx, up_idx] = HPDTile.VERT
                        else:
                            lookup[left_idx, up_idx] = HPDTile.ELBOW_SE

        cls._TBD_LOOKUP = lookup

    @classmethod
    def _get_tbd_tile(cls, left_tile: HPDTile | None, up_tile: HPDTile | None) -> HPDTile:
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
        if HPD._TBD_LOOKUP is None:
            HPD._init_tbd_lookup()

        # Find all TBD positions
        tbd_mask = self._grid == HPDTile.TBD

        # Corner case [0,0]
        if tbd_mask[0, 0]:
            self._grid[0, 0] = HPD._TBD_LOOKUP[8, 8]  # (None, None)

        # First column [1:, 0] - must process sequentially as each row depends on previous
        if rows > 1:
            tbd_rows = np.where(tbd_mask[1:, 0])[0] + 1
            for row in tbd_rows:
                up_tile = int(self._grid[row - 1, 0])
                self._grid[row, 0] = HPD._TBD_LOOKUP[8, up_tile]

        # Process column by column (dependencies require sequential processing)
        for col in range(1, cols):
            # Top row [0, col] - single lookup with (left_tile, None)
            if tbd_mask[0, col]:
                left_tile = int(self._grid[0, col - 1])
                self._grid[0, col] = HPD._TBD_LOOKUP[left_tile, 8]

            # Interior cells [1:, col] - must process sequentially as each row depends on previous
            if rows > 1:
                tbd_rows = np.where(tbd_mask[1:, col])[0] + 1
                for row in tbd_rows:
                    left_tile = int(self._grid[row, col - 1])
                    up_tile = int(self._grid[row - 1, col])
                    self._grid[row, col] = HPD._TBD_LOOKUP[left_tile, up_tile]

    def __len__(self) -> int:
        """Return the size n of the n×n grid"""
        return self.rows

    def __getitem__(self, key) -> HPDTile | np.ndarray:
        """Access grid elements, casting to HPDTile"""
        result = self._grid[key]
        # If it's a single element, convert to HPDTile enum
        if isinstance(result, (int, np.integer)):
            return HPDTile(result)
        # If it's an array, return as-is (caller will need to convert)
        return result

    def _sympyrepr(self, printer=None) -> str:
        """SymPy repr representation"""
        grid_list = self._grid.tolist()
        if printer is None:
            return f"HPD({grid_list})"
        return f"HPD({printer._print(grid_list)})"

    def _sympystr(self, printer=None) -> str:
        """SymPy str representation of the HPD using tile symbols"""
        return printer._print("HPD(\n" + pretty(self) + ")")

    def toggle_bottom_row(self) -> HPD:
        new_id_vector = (*self._id_vector[:-1], 1 - self._id_vector[-1])

        new_grid = self._grid.copy()
        if new_id_vector[-1] == 0:
            weight = np.sum(np.where(self._grid[-1, :] == HPDTile.BLANK, 1, 0))
            new_grid[-1, :] = _classical_bottom_row(weight, self.cols)
        else:
            weight = np.sum(np.where((self._grid[-1, :] == HPDTile.CROSS) | (self._grid[-1, :] == HPDTile.HORIZ), 1, 0))
            new_grid[-1, :] = _bpd_bottom_row(weight, self.cols)
        # print("TOGLLE")
        # _display_grid(new_grid)
        return HPD(new_grid, id_vector=new_id_vector)

    def _pretty(self, printer=None):
        """Pretty printing with row and column labels"""
        from sympy.printing.pretty.stringpict import prettyForm

        if printer is None:
            # Fallback to simple string representation
            return prettyForm(self._sympystr(None))

        # Calculate column widths based on labels (minimum 1 for tile)
        col_widths = []
        for j in range(self.cols):
            label = str(printer._print((~self.perm)[j] if self.perm is not None else "?"))
            col_widths.append(max(1, len(label)))

        # Use maximum column width for all columns
        max_col_width = max(col_widths)

        # Build rows with proper tile extensions
        rows = []
        # perm = self.perm
        # perm_values = [perm[i] for i in range(self.rows)]

        for i in range(self.rows):
            row_parts = []
            for j in range(self.cols):
                tile = self[i, j]
                is_weighty = self.is_weighty_position(i, j)
                # Use context-aware display symbol
                tile_str = tile.get_display_symbol(is_weighty)

                # Symmetric padding to center tile in column using max width
                total_pad = max(max_col_width - 1, 2)
                left_pad = total_pad // 2
                right_pad = total_pad - left_pad

                horiz_str = str(HPDTile.HORIZ)
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
            rows.append(row_str)

        # Calculate max label widths for left and right sides
        max_left_label_width = 0
        max_right_label_width = 0
        for i in range(self.rows):
            row_label = self.row_index_to_label(i)
            row_label_len = len(str(printer._print(row_label)))
            if self._id_vector[i] == 0:
                max_left_label_width = max(max_left_label_width, row_label_len)
            else:
                max_right_label_width = max(max_right_label_width, row_label_len)

        # Now add labels to each row with consistent padding on both sides
        for i in range(len(rows)):
            row_label = str(printer._print(self.row_index_to_label(i)))
            if self._id_vector[i] == 0:
                # Label on the left
                left_part = row_label.ljust(max_left_label_width) + " "
                right_part = " " * (max_right_label_width + 1) if max_right_label_width > 0 else ""
                rows[i] = left_part + rows[i] + right_part
            else:
                # Label on the right
                left_part = " " * (max_left_label_width + 1) if max_left_label_width > 0 else ""
                right_part = " " + row_label.rjust(max_right_label_width)
                rows[i] = left_part + rows[i] + right_part
            # rows.insert(0, f"{row_label} ")  # Add space between rows for readability

        # Add column labels at the top with proper spacing and alignment
        col_labels = []
        # Each column is 1 (tile) + total_pad wide
        total_pad = max(max_col_width - 1, 2)
        column_width = 1 + total_pad

        for j in range(self.cols):
            label = str(printer._print((~self.perm)[j] if self.perm is not None else "?"))
            col_labels.append(label.center(column_width))

        # Add leading and trailing space to align with row labels
        leading_space = " " * (max_left_label_width + 1) if max_left_label_width > 0 else ""
        trailing_space = " " * (max_right_label_width + 1) if max_right_label_width > 0 else ""

        rows.insert(0, leading_space + "".join(col_labels) + trailing_space)

        return prettyForm("\n".join(rows))

    def shiftup(self, shift: int = 1) -> HPD:
        """Shift the HPD up by a given amount."""
        # Create new grid with shifted dimensions
        new_rows = self.rows + shift
        new_cols = self.cols + shift
        new_grid = np.full((new_rows, new_cols), HPDTile.ELBOW_SE, dtype=HPDTile)

        # Shift the grid contents
        for i in range(self.rows):
            for j in range(self.cols):
                new_grid[i + shift, j + shift] = self[i, j]

        # Fill the top-left portion with identity pattern
        for i in range(shift):
            for j in range(shift):
                if i == j:
                    new_grid[i, j] = HPDTile.ELBOW_SE
                elif i < j:
                    new_grid[i, j] = HPDTile.HORIZ
                else:
                    new_grid[i, j] = HPDTile.VERT

        # Connect the identity to shifted content
        for i in range(shift):
            for j in range(shift, new_cols):
                new_grid[i, j] = HPDTile.HORIZ
        for i in range(shift, new_rows):
            for j in range(shift):
                new_grid[i, j] = HPDTile.VERT
        return HPD(new_grid)

    @property
    def perm(self) -> Permutation:
        """
        Compute the permutation associated with this HPD.

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
        # left_col = self._grid[:, 0]
        # right_col = self._grid[:, -1]

        # Check for TBD in bottom row
        # if np.any(bottom_row == HPDTile.TBD):
        #     raise ValueError("Cannot compute permutation with unresolved TBD tiles")

        # Vectorized: find columns with entrance_from_bottom
        # entrance_from_bottom is True for VERT, CROSS, ELBOW_NW
        # left_entrance_mask = (left_col == HPDTile.UHORIZ) | (left_col == HPDTile.WHORIZ) | (left_col == HPDTile.UCROSS) | (left_col == HPDTile.WCROSS) | (left_col == HPDTile.ELBOW_NW) | (left_col == HPDTile.BUMP)
        # right_entrance_mask = (right_col == HPDTile.UHORIZ) | (right_col == HPDTile.WHORIZ) | (right_col == HPDTile.UCROSS) | (right_col == HPDTile.WCROSS) | (right_col == HPDTile.ELBOW_NE)
        # good_rows = (np.where(left_mask)[0] + 1).tolist()

        # good_cols = Permutation.from_partial(good_cols)
        # left_rows = (np.where(left_entrance_mask)[0]).tolist()
        # right_rows = (np.where(right_entrance_mask)[0]).tolist()
        left_spots = [i for i in range(len(self._id_vector)) if self._id_vector[i] == 0]
        left_labels = [self.row_index_to_label(i) for i in left_spots]
        left_rows = list(zip(left_spots, left_labels))
        right_spots = [i for i in range(len(self._id_vector)) if self._id_vector[i] == 1]
        right_labels = [self.row_index_to_label(i) for i in right_spots]
        right_rows = list(zip(right_spots, right_labels))
        top_pop = [0] * self.cols

        def _new_direction(this_row, col, tile, going_right, going_up):
            if tile in (HPDTile.WCROSS, HPDTile.UCROSS):
                if going_up:
                    this_row -= 1
                elif going_right:
                    col += 1
                else:
                    col -= 1
            elif tile in (HPDTile.WHORIZ, HPDTile.UHORIZ):
                if going_right:
                    col += 1
                else:
                    col -= 1
                going_up = False
            elif tile == HPDTile.ELBOW_NW:
                this_row -= 1
                going_up = True
                going_right = False
            elif tile == HPDTile.ELBOW_SE:
                col += 1
                going_up = False
                going_right = True
            elif tile == HPDTile.ELBOW_NE:
                this_row -= 1
                going_right = False
                going_up = True
            elif tile == HPDTile.ELBOW_SW:
                col -= 1
                going_right = False
                going_up = False
            elif tile == HPDTile.VERT:
                this_row -= 1
                going_up = True
                going_right = False
            elif tile == HPDTile.BUMP:
                if going_up:
                    col += 1
                    going_right = True
                    going_up = False
                elif going_right:
                    this_row -= 1
                    going_up = True
                    going_right = False
                else:
                    # print("BADAAD")
                    raise ValueError("Invalid BUMP direction")
            else:
                print(f"Invalid tile {tile} encountered during tracing")
                return None, None, None, None
            return this_row, col, going_right, going_up

        for row, label in left_rows:
            # trace row
            tile = HPDTile(self._grid[row, 0])
            col = 0
            this_row = row
            # print(f"ROW: {row}")
            going_right = True
            going_up = False

            while this_row >= 0:
                # Switch on all tile types
                # print(f"At ({this_row}, {col}) tile {tile}")
                tile = HPDTile(self._grid[this_row, col])
                this_row, col, going_right, going_up = _new_direction(this_row, col, tile, going_right, going_up)
                if this_row is None:
                    self._perm = None
                    return self._perm
            top_pop[col] = label

        for row, label in right_rows:
            # trace row
            col = self.cols - 1
            this_row = row
            going_right = False
            going_up = False
            while this_row >= 0:
                # Switch on all tile types
                tile = HPDTile(self._grid[this_row, col])
                this_row, col, going_right, going_up = _new_direction(this_row, col, tile, going_right, going_up)
                if this_row is None:
                    self._perm = None
                    return self._perm
            top_pop[col] = label
        # print(top_pop)
        # small_perm = Permutation([])
        # Vectorized: Map tiles to their diff values
        # diff = np.ones((nrows, ncols), dtype=int)
        # diff[self._grid == HPDTile.BLANK] = 0
        # diff[self._grid == HPDTile.CROSS] = 2
        # # Create r array with shape (nrows+1, ncols+1)
        # r = np.zeros((nrows + 1, ncols + 1), dtype=int)

        # for i in range(1, nrows + 1):
        #     for j in range(1, ncols + 1):
        #         r[i, j] = r[i - 1, j - 1] + diff[i - 1, j - 1]

        # # Pre-compute all cross positions and their pipes_northeast values
        # cross_positions = np.argwhere(self._grid == HPDTile.CROSS)
        # if len(cross_positions) > 0:
        #     # Sort by column first, then by row descending (for correct swap order)
        #     sort_indices = np.lexsort((-cross_positions[:, 0], cross_positions[:, 1]))
        #     cross_positions = cross_positions[sort_indices]
        #     # Get pipes_northeast for each cross position
        #     pipes_northeast_values = r[cross_positions[:, 0] + 1, cross_positions[:, 1] + 1]
        #     # Apply swaps sequentially
        #     for pipes_northeast in pipes_northeast_values:
        #         small_perm = small_perm.swap(pipes_northeast - 2, pipes_northeast - 1)

        # build_perm = good_cols * small_perm
        print(top_pop)

        self._perm = ~Permutation(top_pop)
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
        Compute the length vector of the permutation represented by this HPD.

        The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
        of weighty tiles in row i. Which tiles are weighty depends on _id_vector[i]:
        - If _id_vector[i] == 1: BLANK tiles are weighty
        - If _id_vector[i] == 0: CROSS and HORIZ tiles are weighty

        Returns:
            Tuple of integers representing the length vector
        """
        # Vectorized approach: different tiles count based on _id_vector
        id_vec = np.array(self._id_vector)[:, np.newaxis]  # Shape (n, 1) for broadcasting

        # Create masks for different tile types
        is_blank = self._grid == HPDTile.BLANK
        is_cross = self._grid == HPDTile.CROSS
        is_horiz = self._grid == HPDTile.HORIZ

        # For id_vector == 1 rows: count BLANK tiles
        # For id_vector == 0 rows: count CROSS and HORIZ tiles
        weighty_mask = np.where(id_vec == 1, is_blank, is_cross | is_horiz)

        # Sum along each row
        weighted_counts = np.sum(weighty_mask, axis=1)
        return tuple(int(weighted_counts[self.row_index_to_label(i) - 1]) for i in range(len(weighted_counts)))

    @classmethod
    def from_asm(cls, asm) -> HPD:
        """
        Create a HPD from an ASM (Alternating Sign Matrix).

        Args:
            asm: n×n array-like of integers (-1, 0, 1)
        Returns:
            HPD object
        """
        corner_sum = np.cumsum(np.cumsum(asm, axis=0), axis=1)
        n = asm.shape[0]
        grid = np.full((n, n), fill_value=HPDTile.TBD, dtype=HPDTile)

        def asm_tile(i, j):
            if corner_sum[i, j] - (corner_sum[i - 1, j - 1] if i > 0 and j > 0 else 0) == 2:
                return HPDTile.CROSS
            if corner_sum[i, j] - (corner_sum[i - 1, j - 1] if i > 0 and j > 0 else 0) == 0:
                return HPDTile.BLANK
            if (corner_sum[i - 1, j] if i > 0 else 0) - (corner_sum[i, j - 1] if j > 0 else 0) == -1:
                return HPDTile.HORIZ
            if (corner_sum[i - 1, j] if i > 0 else 0) - (corner_sum[i, j - 1] if j > 0 else 0) == 1:
                return HPDTile.VERT

            raise ValueError(f"Invalid ASM corner sum at ({i}, {j})")

        asm = np.array(asm, dtype=int)
        grid[asm == 1] = HPDTile.ELBOW_SE
        grid[asm == -1] = HPDTile.ELBOW_NW
        rows, cols = np.where(asm == 0)
        grid[rows, cols] = np.vectorize(asm_tile, otypes=[HPDTile])(rows, cols)
        return cls(grid)

    def to_asm(self):
        to_trans = self.resize(len(self.perm))
        n = max(to_trans.rows, to_trans.cols)
        asm = np.full(shape=(n, n), fill_value=0, dtype=int)
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
    def rothe_bpd(cls, perm: Permutation, num_rows: int | None = None) -> HPD:
        if num_rows is None:
            num_rows = len(perm)
        n = max(num_rows, len(perm))
        grid = np.full((num_rows, n), fill_value=HPDTile.TBD, dtype=HPDTile)

        # Set blanks from diagram
        diagram = perm.diagram
        if diagram:
            diagram_arr = np.array(list(diagram), dtype=int)
            grid[diagram_arr[:, 0] - 1, diagram_arr[:, 1] - 1] = HPDTile.BLANK

        # Vectorize cross detection using graph coordinates
        graph = perm.graph
        if graph and len(graph) > 0:
            graph_arr = np.array(list(graph), dtype=int)

            # Get mask of non-blank positions
            non_blank_mask = grid != HPDTile.BLANK
            non_blank_coords = np.argwhere(non_blank_mask)

            if len(non_blank_coords) > 0:
                # For each non-blank position, check graph conditions
                for i, j in non_blank_coords:
                    # Check if any graph point has row == i+1 and col < j+1
                    has_row_below = np.any((graph_arr[:, 0] == i + 1) & (graph_arr[:, 1] - 1 < j))
                    # Check if any graph point has row < i+1 and col == j+1
                    has_col_left = np.any((graph_arr[:, 0] - 1 < i) & (graph_arr[:, 1] == j + 1))
                    if has_row_below and has_col_left:
                        grid[i, j] = HPDTile.CROSS

        return HPD(grid)

    @property
    def weight(self) -> Tuple[int, ...]:
        """
        Compute the weight of this HPD.

        The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
        of empty squares (0s) in column i.

        Returns:
            Tuple of integers representing the weight
        """
        return tuple(int(np.sum(self[:, j] == HPDTile.BLANK)) for j in range(self.rows))

    @property
    def word(self) -> Tuple[int, ...]:
        """
        Compute a reduced word for the permutation represented by this HPD.

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
        diff[self._grid == HPDTile.BLANK] = 0
        diff[self._grid == HPDTile.CROSS] = 2
        # Create r array with shape (nrows+1, ncols+1)
        r = np.zeros((nrows + 1, ncols + 1), dtype=int)

        for i in range(1, nrows + 1):
            for j in range(1, ncols + 1):
                r[i, j] = r[i - 1, j - 1] + diff[i - 1, j - 1]
        # return r[target_row + 1, target_col + 1]
        for col in range(self.cols):
            for row in range(self.rows - 1, -1, -1):
                if self[row, col] == HPDTile.CROSS:
                    pipes_northeast = r[row + 1, col + 1]  # self.cols] - r[row + 1, col]
                    word.append(pipes_northeast - 1)
        self._word = tuple(word)
        return self._word

    def set_width(self, width):
        """Set the width of the HPD by adding empty columns on the right if needed."""
        # if width < self.cols:
        #     if len(self.perm) > width:
        #         raise ValueError("New width must be at least the length of the permutation")
        #     new_grid = np.full((self.rows, width), fill_value=HPDTile.TBD, dtype=HPDTile)
        #     new_grid[:, :width] = self._grid[:, :width]
        #     bop = HPD(new_grid)
        #     bop.rebuild()
        #     return bop
        # if width == self.cols:
        #     return self
        # new_grid = np.full((self.rows, width), fill_value=HPDTile.TBD, dtype=HPDTile)
        # new_grid[:, : self.cols] = self._grid
        # return HPD(new_grid)
        raise NotImplementedError("set_width is not implemented yet")

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
        if len(self.all_blanks()) != self.perm.inv:
            self._valid = False
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
        """Check equality of two HPDs"""
        if not isinstance(other, HPD):
            return False
        return np.array_equal(self._grid, other._grid)

    def __hash__(self) -> int:
        """Hash for use in sets and dicts"""
        return hash(self._grid.tobytes())

    def copy(self) -> HPD:
        """Create a copy of this HPD"""
        new_bpd = HPD(None, _is_copy=True)
        new_bpd._grid = self._grid.copy()
        new_bpd._column_perm = self._column_perm
        new_bpd._perm = self._perm
        new_bpd._valid = self._valid
        new_bpd._word = self._word
        new_bpd._unzero_cache = self._unzero_cache
        return new_bpd

    @property
    def num_crossings(self) -> int:
        """Total number of crossings in the HPD"""
        return int(np.sum(self._grid == HPDTile.CROSS))

    def right_root_at(self, i: int, j: int) -> int:
        """
        Compute the inversion associated with the crossing at position (i, j).

        The inversion is determined by tracing the pipes through the HPD.

        Args:
            i: Row index of the crossing
            j: Column index of the crossing
        Returns:
            The inversion value as an integer
        """
        if self[i, j] not in (HPDTile.CROSS, HPDTile.BUMP):
            return None
        bpd = self.resize(len(self.perm))
        nrows, ncols = bpd._grid.shape
        # Map tiles to their diff values
        diff = np.ones_like(bpd._grid, dtype=int)
        diff[bpd._grid == HPDTile.BLANK] = 0
        diff[bpd._grid == HPDTile.CROSS] = 2
        # Create r array with shape (nrows+1, ncols+1)
        r = np.zeros((nrows + 1, ncols + 1), dtype=int)

        a, b = 0, 0

        for ai in range(1, nrows + 1):
            for aj in range(1, ncols + 1):
                r[ai, aj] = r[ai - 1, aj - 1] + diff[ai - 1, aj - 1]
        # return r[target_row + 1, target_col + 1]
        for col in range(j, bpd.cols):
            for row in range(bpd.rows - 1, -1, -1):
                if row > i and col == j:
                    continue
                if row == i and col == j:
                    a, b = r[row + 1, col + 1] - 1, r[row + 1, col + 1]
                    continue
                if bpd[row, col] == HPDTile.CROSS:
                    a, b = Permutation.ref_product(r[row + 1, col + 1] - 1).act_root(a, b)  # self.cols] - r[row + 1, col]
        return a, b

    def left_root_at(self, i: int, j: int) -> int:
        """
        Compute the inversion associated with the crossing at position (i, j).

        The inversion is determined by tracing the pipes through the HPD.

        Args:
            i: Row index of the crossing
            j: Column index of the crossing
        Returns:
            The inversion value as an integer
        """
        if self[i, j] not in (HPDTile.CROSS, HPDTile.BUMP):
            return None

        nrows, ncols = self._grid.shape
        # Map tiles to their diff values
        diff = np.ones_like(self._grid, dtype=int)
        diff[self._grid == HPDTile.BLANK] = 0
        diff[self._grid == HPDTile.CROSS] = 2
        # Create r array with shape (nrows+1, ncols+1)
        r = np.zeros((nrows + 1, ncols + 1), dtype=int)

        a, b = 0, 0

        for ai in range(1, nrows + 1):
            for aj in range(1, ncols + 1):
                r[ai, aj] = r[ai - 1, aj - 1] + diff[ai - 1, aj - 1]
        # return r[target_row + 1, target_col + 1]
        for col in range(j, -1, -1):
            for row in range(self.rows):
                if row < i and col == j:
                    continue
                if row == i and col == j:
                    a, b = r[row + 1, col + 1] - 1, r[row + 1, col + 1]
                    continue
                if self[row, col] == HPDTile.CROSS:
                    a, b = Permutation.ref_product(r[row + 1, col + 1] - 1).act_root(a, b)  # self.cols] - r[row + 1, col]
        return a, b
        # word.append(pipes_northeast - 1)

    # def inversion_at_bump(self, i: int, j: int) -> int:
    #     """
    #     Compute the inversion associated with the crossing at position (i, j).

    #     The inversion is determined by tracing the pipes through the HPD.

    #     Args:
    #         i: Row index of the crossing
    #         j: Column index of the crossing
    #     Returns:
    #         The inversion value as an integer
    #     """
    #     if self[i, j] != HPDTile.BUMP:
    #         return None

    #     # up and right
    #     up_r = i - 1
    #     up_c = j

    #     while up_r >= 0 and up_c < self.cols:
    #         tile = self[up_r, up_c]
    #         if tile.feeds_up:
    #             up_r -= 1
    #         else:
    #             up_c += 1

    #     right_r = i
    #     right_c = j + 1

    #     while right_r >= 0 and right_c < self.cols:
    #         tile = self[right_r, right_c]
    #         if tile.feeds_right:
    #             right_c += 1
    #         else:
    #             right_r -= 1
    #     assert up_r < right_r
    #     return up_r + 1, right_r + 1

    def trace_pipe(self, i: int, j: int, direction: str | None = None) -> int | None:
        if self[i, j] == HPDTile.ELBOW_NW:
            if direction == "left":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self[i, j] == HPDTile.ELBOW_SE:
            if direction == "down":
                return None
            if i == self.rows - 1:
                return self._column_perm[j]
            return self.trace_pipe(i + 1, j, direction="down")
        if self[i, j] == HPDTile.HORIZ:
            if direction == "down":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self[i, j] == HPDTile.VERT:
            if direction == "left":
                return None
            if i == self.rows - 1:
                return self._column_perm[j]
            return self.trace_pipe(i + 1, j, direction="down")
        if self[i, j] == HPDTile.CROSS:
            if direction == "down":
                if i == self.rows - 1:
                    return self._column_perm[j]
                return self.trace_pipe(i + 1, j, direction="down")
            if direction == "left":
                return self.trace_pipe(i, j - 1, direction="left")
            raise ValueError("Must specify direction when tracing through a crossing")
        if self[i, j] == HPDTile.BUMP:
            if direction == "down":
                return self.trace_pipe(i, j - 1, direction="left")
            if direction == "left":
                if i == self.rows - 1:
                    return self._column_perm[j]
                return self.trace_pipe(i + 1, j, direction="down")
            raise ValueError("Must specify direction when tracing through a crossing")
        raise ValueError(f"Invalid tile for tracing pipe at ({i}, {j}): {self[i, j]}\n{self=}")

    def all_se_elbows(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(HPDTile.ELBOW_SE)

    def all_nw_elbows(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(HPDTile.ELBOW_NW)

    def all_blanks(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(HPDTile.BLANK)

    def all_crossings(self) -> set[tuple[int, int]]:
        return self.all_tiles_of_type(HPDTile.CROSS)

    def all_tiles_of_type(self, tile_type: HPDTile) -> set[tuple[int, int]]:
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
                    if (i, j) != (ri, rj) and (self[i, j] == HPDTile.ELBOW_SE or self[i, j] == HPDTile.ELBOW_NW or self[i, j] == HPDTile.BUMP):  # if not NW-corner, check if elbow
                        LEGAL = False
                        break
                if LEGAL:
                    droop_moves.add(((ri, rj), (bi, bj)))
        return droop_moves

    def min_droop_moves(self) -> set[tuple[tuple[int, int], tuple[int, int]]]:
        droop_moves = set()
        for a, b in self.all_tiles_of_type((HPDTile.ELBOW_SE, HPDTile.BUMP)):
            # Find min x coordinate in column b where tile is not CROSS and x > a
            non_cross_positions = np.argwhere((self._grid[:, b] != HPDTile.CROSS) & (np.arange(self.rows) > a))
            if len(non_cross_positions) > 0:
                x = non_cross_positions[0][0]
            else:
                continue
            non_cross_positions = np.argwhere((self._grid[a, :] != HPDTile.CROSS) & (np.arange(self.cols) > b))
            # y = np.min(np.argwhere(self._grid[a, b + 1 :] != HPDTile.CROSS), initial=ARBITRARY_LARGE) + b + 1
            if len(non_cross_positions) > 0:
                y = non_cross_positions[0][0]
            else:
                continue
            droop_moves.add(((a, b), (x, y)))
        return droop_moves

    def do_min_droop_move(self, move: tuple[tuple[int, int], tuple[int, int]]) -> HPD:
        D = self.copy()
        (ri, rj) = move[0]
        (bi, bj) = move[1]

        # if bj >= self.cols or bi >= self.rows:
        #     D = self.resize(max(bi + 2, bj + 2))
        orig_self = self.copy()
        D._grid[ri, rj] = HPDTile.BLANK if orig_self[ri, rj] == HPDTile.ELBOW_SE else HPDTile.ELBOW_NW  # NW-corner
        D._grid[bi, bj] = HPDTile.ELBOW_NW if orig_self[bi, bj] == HPDTile.BLANK else HPDTile.BUMP  # SE-corner
        D._grid[bi, rj] = HPDTile.ELBOW_SE  # SW-corner
        D._grid[ri, bj] = HPDTile.ELBOW_SE  # NE-corner

        # top and bottom
        for j in range(rj + 1, bj):
            if orig_self[ri, j] == HPDTile.HORIZ:
                D._grid[ri, j] = HPDTile.BLANK
            else:  # self[ri, j] == HPDTile.CROSS
                D._grid[ri, j] = HPDTile.VERT
            if orig_self[bi, j] == HPDTile.BLANK:
                D._grid[bi, j] = HPDTile.HORIZ
            else:  # self[bi, j] == HPDTile.VERT
                D._grid[bi, j] = HPDTile.CROSS
        # left and right
        for i in range(ri + 1, bi):
            if orig_self[i, rj] == HPDTile.VERT:
                D._grid[i, rj] = HPDTile.BLANK
            else:  # self[i, rj] == HPDTile.CROSS
                D._grid[i, rj] = HPDTile.HORIZ
            if orig_self[i, bj] == HPDTile.BLANK:
                D._grid[i, bj] = HPDTile.VERT
            else:
                D._grid[i, bj] = HPDTile.CROSS
        return D

    def do_droop_move(self, move: tuple[tuple[int, int], tuple[int, int]]) -> HPD:
        D = self.copy()
        (ri, rj) = move[0]
        (bi, bj) = move[1]

        D._grid[ri, rj] = HPDTile.BLANK  # NW-corner
        D._grid[bi, bj] = HPDTile.ELBOW_NW  # SE-corner
        D._grid[bi, rj] = HPDTile.ELBOW_SE  # SW-corner
        D._grid[ri, bj] = HPDTile.ELBOW_SE  # NE-corner

        # top and bottom
        for j in range(rj + 1, bj):
            if self[ri, j] == HPDTile.HORIZ:
                D._grid[ri, j] = HPDTile.BLANK
            else:  # self[ri, j] == HPDTile.CROSS
                D._grid[ri, j] = HPDTile.VERT
            if self[bi, j] == HPDTile.BLANK:
                D._grid[bi, j] = HPDTile.HORIZ
            else:  # self[bi, j] == HPDTile.VERT
                D._grid[bi, j] = HPDTile.CROSS
        # left and right
        for i in range(ri + 1, bi):
            if self[i, rj] == HPDTile.VERT:
                D._grid[i, rj] = HPDTile.BLANK
            else:  # self[i, rj] == HPDTile.CROSS
                D._grid[i, rj] = HPDTile.HORIZ
            if self[i, bj] == HPDTile.BLANK:
                D._grid[i, bj] = HPDTile.VERT
            else:  # self[i, bj] == HPDTile.HORIZ
                D._grid[i, bj] = HPDTile.CROSS
        D.rebuild()
        return D

    def huang_bump(self, a, b):
        # perm_inverse = ~self.perm
        working_bpd = self.copy()  # resize(len(self.perm))
        the_cross_list = [(i, j) for (i, j) in working_bpd.all_crossings() if working_bpd.right_root_at(i, j) == (a, b)]
        if len(the_cross_list) == 0:
            raise ValueError(f"No crossing found for inversion ({a}, {b})")
        x, y = the_cross_list[0]
        working_bpd._grid[x, y] = HPDTile.BUMP
        working_bpd.rebuild()
        return working_bpd._monk_iterate(x, y, x + 1, self.perm.swap(a - 1, b - 1)[x])

    @classmethod
    def random_bpd(cls, perm, num_rows):
        import random

        return random.choice(list(cls.all_bpds(perm, num_rows)))

    def monk_insert(self, row):
        """RETURNS NORMALIZED"""
        # find easternmost tile
        working_bpd = self.copy()  # resize(max(row, len(self.perm), self.rows))
        # if len(working_bpd) < row:
        #     working_bpd = working_bpd.resize(row)
        alpha = row - 1
        col_coord = np.max(np.argwhere(working_bpd._grid[alpha, :] == HPDTile.ELBOW_SE))
        x, y = alpha, col_coord
        return working_bpd._monk_iterate(x, y, row, self.perm[row - 1]).resize(max(row, self.rows))

    def _monk_iterate(self, x, y, row, row_val) -> HPD:
        x_iter, y_iter = x, y
        working_bpd = self.copy()
        while True:
            min_droops = working_bpd.min_droop_moves()
            try:
                the_move = next(((a, b), (c, d)) for ((a, b), (c, d)) in min_droops if (a == x_iter) and (b == y_iter))
            except StopIteration:
                # print(f"No droop move found for position ({x_iter}, {y_iter}) in HPD:\n{working_bpd}")
                # print(f"{min_droops=}")
                raise
            working_bpd = working_bpd.do_min_droop_move(the_move)
            if working_bpd[the_move[1][0], the_move[1][1]] == HPDTile.ELBOW_NW:
                spots = [(a, b) for (a, b) in working_bpd.all_se_elbows() if working_bpd.trace_pipe(a, b) == row_val and a == the_move[1][0]]
                if len(spots) == 0:
                    raise ValueError(f"No SE elbow found for pipe {row} after droop move in monk iteration.")
                assert len(spots) == 1
                x_iter, y_iter = spots[0]
                continue
            z, w = the_move[1]
            assert working_bpd[z, w] == HPDTile.BUMP
            pipe1, pipe2 = working_bpd.right_root_at(z, w)
            any_cross = [(zp, wp) for (zp, wp) in working_bpd.all_crossings() if set(working_bpd.left_root_at(zp, wp)) == {pipe1, pipe2} and zp != z and wp != w]
            if any_cross:
                z_prime, w_prime = any_cross[0]
                working_bpd._grid[z, w], working_bpd._grid[z_prime, w_prime] = working_bpd._grid[z_prime, w_prime], working_bpd._grid[z, w]
                working_bpd.rebuild()
                x_iter, y_iter = z_prime, w_prime
                continue
            working_bpd._grid[z, w] = HPDTile.CROSS
            break
        working_bpd.rebuild()
        return working_bpd

    def normalize(self) -> HPD:
        raise NotImplementedError("normalize method is not implemented yet")
        # if len(self) == len(self.perm):
        #     if len(self.perm) == self.cols:
        #         return self.copy()
        #     return self.set_width(len(self.perm))
        # snap_size = len(self.perm)
        # new_grid = self._grid.copy()
        # if snap_size < self.cols:
        #     if snap_size < len(self):
        #         new_grid = new_grid[:snap_size, :snap_size]
        #     else:
        #         new_grid = new_grid[:, :snap_size]
        # if snap_size > new_grid.shape[0] or snap_size > new_grid.shape[1]:
        #     new_grid = np.pad(new_grid, ((0, snap_size - new_grid.shape[0]), (0, max(0, snap_size - new_grid.shape[1]))), constant_values=HPDTile.TBD)
        # bottom_portion = HPD.rothe_bpd(self.perm.min_coset_rep(*(list(range(self.rows)) + list(range(self.rows + 1, snap_size)))), snap_size)
        # new_grid[self.rows :, :] = bottom_portion._grid[self.rows :, :]
        # _invalidate_grid(new_grid)
        # new_bpd = HPD(new_grid)
        # assert new_bpd.rows == new_bpd.cols == snap_size
        # return new_bpd

    def pop_op(self) -> tuple[HPD, tuple[int, int]]:
        # --- STEP 0 --- #
        # if len(self) < len(self.perm):
        #     # D = self.normalize()
        #     D = self.resize(len(self.perm))
        # else:
        D = self.resize(max(self.rows, len(self.perm)))

        assert D.is_reduced
        # check if D has a blank tile (i.e., the coxeter length of D.w is zero)
        # if self.perm.inv == 0:
        #     raise ValueError("Cannot perform pop_op on a HPD with zero Coxeter length")

        # find the first row r with a blank tile - vectorized
        blank_positions = np.argwhere(D._grid == HPDTile.BLANK)
        # if len(blank_positions) == 0:
        #     raise ValueError("Cannot perform pop_op on a HPD with no blanks")
        r = blank_positions[0, 0]

        # initialize the mark X at the blank tile (x,y)
        x = r
        y = np.max(blank_positions[blank_positions[:, 0] == r, 1])

        while True:
            # --- STEP 1 --- #
            # move the mark to the rightmost blank tile in the block - vectorized
            row_blanks = D._grid[x, y : D.cols] == HPDTile.BLANK
            if np.any(row_blanks):
                y = y + np.where(~row_blanks)[0][0] - 1 if np.any(~row_blanks) else D.cols - 1

            # --- STEP 2 --- #
            # find first j-elbow (x_,y+1) with x_>x - vectorized search
            if y + 1 >= D.cols:
                break
            elbow_positions = np.where(D._grid[x + 1 : D.rows, y + 1] == HPDTile.ELBOW_NW)[0]
            x_ = elbow_positions[0] + x + 1 if len(elbow_positions) > 0 else 0

            if x_ == 0:  # p==y+1
                break

            # Vectorized tile transformation for z in range(x+1, x_)
            z_range = slice(x + 1, x_)
            col_y_tiles = D._grid[z_range, y]

            # Create masks for each tile type
            blank_mask = col_y_tiles == HPDTile.BLANK
            horiz_mask = col_y_tiles == HPDTile.HORIZ
            elbow_se_mask = col_y_tiles == HPDTile.ELBOW_SE
            elbow_nw_mask = col_y_tiles == HPDTile.ELBOW_NW

            # Apply transformations
            D._grid[z_range, y][blank_mask] = HPDTile.VERT
            D._grid[z_range, y + 1][blank_mask] = HPDTile.BLANK
            D._grid[z_range, y][horiz_mask] = HPDTile.CROSS
            D._grid[z_range, y + 1][horiz_mask] = HPDTile.HORIZ
            D._grid[z_range, y][elbow_se_mask] = HPDTile.VERT
            D._grid[z_range, y + 1][elbow_se_mask] = HPDTile.ELBOW_SE
            D._grid[z_range, y][elbow_nw_mask] = HPDTile.CROSS
            D._grid[z_range, y + 1][elbow_nw_mask] = HPDTile.ELBOW_NW

            D._grid[x, y] = HPDTile.ELBOW_SE  # NW-corner
            D._grid[x_, y + 1] = HPDTile.BLANK  # SE-corner
            if D[x_, y] == HPDTile.ELBOW_SE:  # SW-corner
                D._grid[x_, y] = HPDTile.VERT
            else:  # D._grid[x_,y] == HPDTile.HORIZ
                D._grid[x_, y] = HPDTile.ELBOW_NW
            if D[x, y + 1] == HPDTile.ELBOW_SE:  # NE-corner
                D._grid[x, y + 1] = HPDTile.HORIZ
            else:  # D._grid[x,y+1] == HPDTile.VERT
                D._grid[x, y + 1] = HPDTile.ELBOW_NW

            # move the mark X to the SE-corner of U
            x = x_
            y = y + 1

        # --- STEP 3 --- #
        a = y  # where (x,y) is the final position of the mark X

        # Find x_ - vectorized search from bottom
        if y < D.cols:
            elbow_positions = np.where(D._grid[x + 1 : D.rows, y] == HPDTile.ELBOW_SE)[0]
            x_ = elbow_positions[-1] + x + 1 if len(elbow_positions) > 0 else 0
        else:
            x_ = 0

        # Vectorized tile transformation for z in range(x+1, x_)
        if x_ > x + 1:
            z_range = slice(x + 1, x_)
            col_y_tiles = D._grid[z_range, y]

            # Create masks for each tile type
            blank_mask = col_y_tiles == HPDTile.BLANK
            horiz_mask = col_y_tiles == HPDTile.HORIZ
            elbow_se_mask = col_y_tiles == HPDTile.ELBOW_SE
            elbow_nw_mask = col_y_tiles == HPDTile.ELBOW_NW

            # Apply transformations
            D._grid[z_range, y][blank_mask] = HPDTile.VERT
            D._grid[z_range, y + 1][blank_mask] = HPDTile.BLANK
            D._grid[z_range, y][horiz_mask] = HPDTile.CROSS
            D._grid[z_range, y + 1][horiz_mask] = HPDTile.HORIZ
            D._grid[z_range, y][elbow_se_mask] = HPDTile.VERT
            D._grid[z_range, y + 1][elbow_se_mask] = HPDTile.ELBOW_SE
            D._grid[z_range, y][elbow_nw_mask] = HPDTile.CROSS
            D._grid[z_range, y + 1][elbow_nw_mask] = HPDTile.ELBOW_NW

        if x_ > 0:
            D._grid[x, y] = HPDTile.ELBOW_SE  # NW-corner
            D._grid[x_, y + 1] = HPDTile.ELBOW_SE  # SE-corner
            D._grid[x_, y] = HPDTile.VERT  # SW-corner
            if D[x, y + 1] == HPDTile.ELBOW_SE:  # NE-corner
                D._grid[x, y + 1] = HPDTile.HORIZ
            else:  # D._grid[x,y+1] == HPDTile.VERT
                D._grid[x, y + 1] = HPDTile.ELBOW_NW

        D = D.resize(self.rows)
        if self.DEBUG:
            assert D.perm.inv == self.perm.inv - 1, f"Resulting HPD inversion count incorrect after pop_op: {D.perm.inv} vs {self.perm.inv - 1} \n{D}\n{self=}"
        return D, (int(a + 1), int(r + 1))

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

    def resize(self, new_num_rows: int) -> HPD:
        self._perm = None  # Invalidate cached permutation
        if new_num_rows > self.rows:
            new_bpd = self.copy()
            # if new_num_rows > self.cols:
            #     new_bpd.set_width(new_num_rows)
            new_grid = np.pad(new_bpd._grid, ((0, new_num_rows - new_bpd.rows), (0, max(0, max(new_num_rows, len(self.perm)) - new_bpd.cols))), constant_values=HPDTile.TBD)
            # elif new_bpd.rows > new_num_rows:
            #     new_bpd._grid = new_bpd._grid[:new_num_rows, : max(new_num_rows, len(self.perm))]
            # bottom_portion = HPD.rothe_bpd(self.perm.min_coset_rep(*(list(range(self.rows)) + list(range(self.rows + 1, max(len(self.perm),new_num_rows)))))
            # new_grid[self.rows :, :] = bottom_portion._grid[self.rows :, :]
            # for r in range(self.rows):
            #     current_col = self.cols
            #     while current_col < new_grid.shape[1]:
            #         if new_grid[r, current_col] == HPDTile.TBD:
            #             new_grid[r, current_col] = HPDTile.HORIZ
            #             current_col += 1
            # _invalidate_grid(new_grid)
            bottom_portion = HPD.rothe_bpd(self.perm, new_num_rows)

            new_grid[self.rows :, : bottom_portion.cols] = bottom_portion._grid[self.rows :, :]
            _invalidate_grid(new_grid)
            new_bpd = HPD(new_grid)
            if new_bpd.cols > max(new_num_rows, len(new_bpd.perm)):
                new_bpd._grid = new_bpd._grid[:, : max(new_num_rows, len(new_bpd.perm))]
                new_bpd.rebuild()
            if self.DEBUG:
                assert new_bpd.is_valid, f"Resulting HPD is not valid after increasing size, \n{pretty(self)} {new_num_rows} {new_bpd!r}"
            return new_bpd
        if new_num_rows < self.rows:
            new_bpd = HPD(self._grid[:new_num_rows, : max(len(self.perm), new_num_rows)])
            new_bpd.rebuild()
            if new_bpd.cols > max(new_num_rows, len(new_bpd.perm)):
                new_bpd._grid = new_bpd._grid[:, : max(new_num_rows, len(new_bpd.perm))]
                new_bpd.rebuild()
            return new_bpd
        new_bpd = self.copy()
        if new_bpd.cols > max(new_num_rows, len(new_bpd.perm)):
            new_bpd._grid = new_bpd._grid[:, : max(new_num_rows, len(new_bpd.perm))]
            new_bpd.rebuild()
        return new_bpd

    def combine(self, other, shift=None) -> HPD:
        """Shift the HPD up by adding empty rows at the bottom."""
        return HPD.from_rc_graph(RCGraph([*self.to_rc_graph()[:shift], *other.to_rc_graph().shiftup(shift)]))

    # def complete_partial(self, partition = None) -> HPD:
    #     """Complete a partial HPD to a full HPD."""
    #     if partition is None:
    #         partition = [self.cols] * self.rows
    #     if len(partition) < len(self.perm):
    #         partition = partition + [0] * (len(self.perm) - len(partition))
    #     new_grid = self._grid.copy()
    #     new_grid = np.pad(new_grid, ((0, self.rows - len(partition)), (0, max(0, len(self.perm) - self.cols))), constant_values=HPDTile.TBD)

    # def left_row_act(self, p: int):
    #     new_grid = self._grid.copy()
    #     new_grid = np.pad(new_grid, ((1, 0), (0,1)), constant_values=HPDTile.TBD)

    def product(self, other: HPD) -> dict[HPD, int]:
        """Compute the product of this HPD with another."""
        from schubmult.utils.perm_utils import add_perm_dict

        if len(other) == 0:
            return {self: 1}
        self_perm = self.perm
        other_shift = other.shiftup(len(self))
        if self_perm.inv == 0:
            # return {HPD.rothe_bpd(Permutation([]), len(self) + len(other)).inverse_pop_op(*other_reduced_compatible).resize(len(self) + len(other)): 1}
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

    def prod_with_bpd(self, other: HPD) -> HPD:
        """Deprecated: Use product() instead. Returns the single HPD from product dictionary."""
        return self.product(other)

    def inverse_pop_op(self, *interlaced_rc) -> HPD:
        if self.rows < len(self.perm):
            D = self.resize(len(self.perm))
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
            elbow_positions = np.where(D._grid[: D.rows, a] == HPDTile.ELBOW_SE)[0]
            if len(elbow_positions) == 0:
                raise ValueError("No elbow found in specified column for inverse pop operation when inserting ")
            x_ = elbow_positions[-1]  # Last (bottom-most) elbow
            D._grid[x_, a] = HPDTile.CROSS
            y = a - 1
            # find x in column y, it will be an SE elbow - vectorized search
            if y >= 0:
                elbow_positions = np.where(D._grid[:x_, y] == HPDTile.ELBOW_SE)[0]
                if len(elbow_positions) == 0:
                    raise ValueError("No elbow found in specified column for inverse pop operation")
                x = elbow_positions[-1]  # Last (bottom-most) elbow before x_
            else:
                x = -1
                raise ValueError("No elbow found in specified column for inverse pop operation")
            D._grid[x, y] = HPDTile.BLANK

            # Vectorized tile transformation for z in range(x+1, x_)
            if x_ > x + 1:
                z_range = slice(x + 1, x_)
                col_y = D._grid[z_range, y]
                col_y1 = D._grid[z_range, y + 1]

                # Create masks for tile patterns
                mask1 = (col_y == HPDTile.VERT) & (col_y1 == HPDTile.BLANK)
                mask2 = (col_y == HPDTile.CROSS) & (col_y1 == HPDTile.ELBOW_NW)
                mask3 = (col_y == HPDTile.CROSS) & (col_y1 == HPDTile.HORIZ)
                mask4 = (col_y == HPDTile.VERT) & (col_y1 == HPDTile.ELBOW_SE)

                # Apply transformations
                D._grid[z_range, y][mask1] = HPDTile.BLANK
                D._grid[z_range, y + 1][mask1] = HPDTile.VERT
                D._grid[z_range, y][mask2] = HPDTile.TBD
                D._grid[z_range, y + 1][mask2] = HPDTile.TBD
                D._grid[z_range, y][mask3] = HPDTile.TBD
                D._grid[z_range, y + 1][mask3] = HPDTile.CROSS
                D._grid[z_range, y][mask4] = HPDTile.ELBOW_SE
                D._grid[z_range, y + 1][mask4] = HPDTile.CROSS

            D.rebuild()

            while True:
                if x == r - 1:
                    break

                # find x_
                x_ = x
                y = y - 1
                if self.DEBUG:
                    assert D[x_, y + 1] == HPDTile.BLANK, "Expected NW elbow during inverse pop operation"
                # Vectorized search: find first non-blank from left in row
                if y >= 0:
                    row_non_blanks = D._grid[x_, : y + 1] != HPDTile.BLANK
                    non_blank_positions = np.where(row_non_blanks)[0]
                    if len(non_blank_positions) > 0:
                        y = non_blank_positions[-1]  # Rightmost non-blank
                    else:
                        y = -1
                else:
                    y = -1
                D._grid[x_, y + 1] = HPDTile.ELBOW_NW

                # find x at SE elbow - vectorized search
                if y >= 0:
                    elbow_positions = np.where(D._grid[:x_, y] == HPDTile.ELBOW_SE)[0]
                    if len(elbow_positions) == 0:
                        raise ValueError("No elbow found in specified column for inverse pop operation")
                    x = elbow_positions[-1]
                else:
                    raise ValueError("No elbow found in specified column for inverse pop operation")

                # [x, y] becomes
                D._grid[x, y] = HPDTile.BLANK

                # Vectorized tile transformation for z in range(x+1, x_)
                if x_ > x + 1:
                    z_range = slice(x + 1, x_)
                    col_y = D._grid[z_range, y]
                    col_y1 = D._grid[z_range, y + 1]

                    # Create masks for tile patterns
                    mask1 = (col_y == HPDTile.VERT) & (col_y1 == HPDTile.BLANK)
                    mask2 = (col_y == HPDTile.CROSS) & (col_y1 == HPDTile.ELBOW_NW)
                    mask3 = (col_y == HPDTile.CROSS) & (col_y1 == HPDTile.HORIZ)
                    mask4 = (col_y == HPDTile.VERT) & (col_y1 == HPDTile.ELBOW_SE)

                    # Apply transformations
                    D._grid[z_range, y][mask1] = HPDTile.BLANK
                    D._grid[z_range, y + 1][mask1] = HPDTile.VERT
                    D._grid[z_range, y][mask2] = HPDTile.TBD
                    D._grid[z_range, y + 1][mask2] = HPDTile.TBD
                    D._grid[z_range, y][mask3] = HPDTile.TBD
                    D._grid[z_range, y + 1][mask3] = HPDTile.CROSS
                    D._grid[z_range, y][mask4] = HPDTile.ELBOW_SE
                    D._grid[z_range, y + 1][mask4] = HPDTile.CROSS

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
            raise ValueError("HPD is not valid, cannot compute reduced compatible sequence")
        if not self.is_reduced:
            raise ValueError("HPD is not reduced, cannot compute reduced compatible sequence")
        work_bpd = self
        ret = []
        # last_inv = work_bpd.perm.inv
        while work_bpd.perm.inv > 0:
            work_bpd, (reflection, row) = work_bpd.pop_op()
            # assert work_bpd.perm.inv == last_inv - 1, "Invariant violated in as_reduced_compatible"
            # last_inv = work_bpd.perm.inv
            ret.append((reflection, row))
        ret.reverse()
        if self.DEBUG:
            assert len(ret) == self.perm.inv
        return tuple(ret)

    def rebuild(self) -> None:
        """Rebuild the HPD to resolve any TBD tiles"""
        # Keep only BLANK and CROSS tiles, wipe out everything else to TBD
        _invalidate_grid(self._grid)
        self._perm = None
        self._valid = None
        self._word = None
        self._unzero_cache = None
        # self.build()

    def zero_out_last_row(self) -> HPD:
        return self.resize(self.rows - 1)

    def set_tile(self, i: int, j: int, tile_type: HPDTile) -> None:
        new_bpd = self.copy()
        new_bpd._grid[i, j] = tile_type
        new_bpd.rebuild()
        return new_bpd

    def right_zero_act(self) -> set[HPD]:
        # # find crosses, untransition them
        if self._unzero_cache is not None:
            return self._unzero_cache
        if not self.is_valid:
            return set()
        resized = self.resize(len(self.perm) + 1)
        results = set()

        crossings = np.array([(i, j) for (i, j) in resized.all_crossings() if i == self.rows], dtype=int)

        n_crossings = len(crossings)
        base_grid = resized._grid.copy()
        for mask in range(1 << n_crossings):  # 2^n combinations
            # Start from base grid instead of full HPD copy
            working_grid = base_grid.copy()
            # Apply mask: set to TBD where mask bit is 1
            for idx in range(n_crossings):
                if mask & (1 << idx):
                    i, j = crossings[idx]
                    working_grid[i, j] = HPDTile.TBD
            working_grid = working_grid[: self.rows + 1, :]
            _invalidate_grid(working_grid)
            new_bpd = HPD(working_grid)

            if new_bpd.is_valid and new_bpd.is_reduced:
                if new_bpd.rows != self.rows + 1 or new_bpd.cols != max(self.rows + 1, len(new_bpd.perm)):
                    new_bpd = new_bpd.resize(self.rows + 1)
                results.add(new_bpd)
        # results = set()

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
        #     new_grid = np.pad(new_grid, ((0, row_pad), (0, col_pad)), constant_values=HPDTile.TBD)
        #     moving_perm =

        self._unzero_cache = results
        return results

    def set_tiles(self, a, b, value: HPDTile) -> None:
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
    #     # Start from base grid instead of full HPD copy
    #     working_grid = base_grid.copy()
    #     # Apply mask: set to TBD where mask bit is 1
    #     for idx in range(n_crossings):
    #         if mask & (1 << idx):
    #             i, j = crossings[idx]
    #             working_grid[i, j] = HPDTile.TBD
    #     _invalidate_grid(working_grid)
    #     new_bpd = HPD(working_grid).snap_width()

    #     if new_bpd.is_valid and new_bpd.is_reduced:
    #         results.add(new_bpd)

    # return results

    def snap_width(self) -> HPD:
        """Snap the width of the HPD to the length of its permutation."""
        perm_length = len(self.perm)
        if self.cols == perm_length or (self.rows > perm_length and self.cols == self.rows):
            return self
        new_bpd = self.copy()
        new_bpd._grid = new_bpd._grid[:, :perm_length]
        # _invalidate_grid(new_bpd._grid)
        # return self.set_width(max(self.rows, perm_length))
        return new_bpd

    def polyvalue(self, x: Sequence[Expr], y: Sequence[Expr] | None = None, **_kwargs) -> Expr:
        """
        Compute the Schubert polynomial value for this HPD.

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
        Convert this HPD to an RC-graph representation.

        Returns:
            RCGraph object (if available in the module)
        """

        rows = [[] for _ in range(self.rows)]
        rc = self.as_reduced_compatible()
        for reflection, row in reversed(rc):
            rows[row - 1] = rows[row - 1] + [reflection]

        return RCGraph([tuple(r) for r in rows])

"""
Bumpless Pipe Dreams (BPD) module

A bumpless pipe dream is a pipe dream diagram with no elbows - only crossings and empty boxes.
Represented as an n x n grid where:
- 1 = crossing (pipes cross)
- 0 = empty box (pipes go straight through)

For general pipe dreams, there are 6 possible tile types.
"""

from enum import IntEnum
from functools import cached_property
from typing import Tuple

import numpy as np

from schubmult.schub_lib.perm_lib import Permutation


class TileType(IntEnum):
    """
    Enumeration of the 6 possible tile types in a pipe dream.

    Each tile represents how two pipes (horizontal and vertical) interact in a square.
    """

    TBD = 0  # Placeholder for uninitialized tile
    EMPTY = 1  # Both pipes go straight (no crossing, no elbow)
    CROSS = 2  # Pipes cross each other
    HORIZ = 3
    ELBOW_NW = 4  # Elbow: bottom-right to top-left (╯)
    ELBOW_SE = 5  # Elbow: top-left to bottom-right (╮)
    VERT = 6
    RC_EMPTY = 7  # Special empty tile for RC-graphs

    def __str__(self):
        symbols = {TileType.EMPTY: "·", TileType.CROSS: "┼", TileType.ELBOW_NW: "╯", TileType.ELBOW_SE: "╭", TileType.HORIZ: "─", TileType.VERT: "│", TileType.RC_EMPTY: "o", TileType.TBD: "?"}
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
        return self == TileType.EMPTY

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


class BPD:
    """
    Bumpless Pipe Dream representation.

    A bumpless pipe dream is an n×n grid where:
    - TileType.CROSS (1) represents a crossing
    - TileType.EMPTY (0) represents an empty box (pipes go straight)
    - For general pipe dreams, can use TileType.ELBOW_* (2-5) for elbows

    Each BPD corresponds to a permutation and has an associated weight.
    """

    def __init__(self, grid):
        """
        Initialize a BPD from a grid.

        Args:
            grid: n×n array-like of TileType values, integers 0-5, or list of lists
        """
        self.grid = np.array(grid, dtype=TileType)

        # Validate grid is square
        if len(self.grid.shape) != 2 or self.grid.shape[0] != self.grid.shape[1]:
            raise ValueError("BPD grid must be square (n×n)")

        self.n = self.grid.shape[0]

        self.build()

    # def _get_tbd_tile(self, left_tile, right_tile, down_tile, up_tile) -> TileType:
    #     # bottom left, vert or elbow
    #     if down_tile is None:
    #         if left_tile is None:
    #             if right_tile.entrance_from_left:
    #                 if up_tile.entrance_from_bottom:
    #                     raise ValueError("Cannot have both right and up tiles not feeding into bottom-left tile")
    #                 return TileType.ELBOW_SE
    #             if up_tile.entrance_from_bottom:
    #                 return TileType.VERT
    #             raise ValueError("At least one of right or up tiles must feed into bottom-left tile")
    #         if left_tile.feeds_right:
    #             if not up_tile.entrance_from_bottom:
    #                 raise ValueError("Cannot have both left and up tiles not feeding into bottom-left tile")
    #             return TileType.CROSS
    #         if right_tile is None or right_tile.entrance_from_left:
    #             if up_tile.entrance_from_bottom:
    #                 raise ValueError("Cannot have both left and up tiles not feeding into bottom-left tile")
    #             return TileType.ELBOW_SE
    #         return TileType.VERT
    #     if left_tile is None:
    #         # must have bottom feed
    #         if up_tile is None or not up_tile.entrance_from_bottom:
    #             return TileType.ELBOW_SE
    #         if right_tile is not None and right_tile.entrance_from_left:
    #             raise ValueError("Cannot have both right and up tiles feeding into bottom-left tile")
    #         return TileType.VERT
    #     if not left_tile.feeds_right and not down_tile.feeds_up:
    #         raise ValueError("Down or left must feed into tile that is not on the bottom or left")
    #     if up_tile is None:
    #         # something must feed in
    #         if right_tile is not None and not right_tile.entrance_from_left:
    #             raise ValueError("Upper TBD tile must exit right")
    #         if left_tile.feeds_right and down_tile.feeds_up:
    #             raise ValueError("Cannot have both left and down tiles feeding into upper tile")
    #         if left_tile.feeds_right:
    #             return TileType.HORIZ
    #         return TileType.ELBOW_SE
    #     if (left_tile.feeds_right and down_tile.feeds_up) and (not up_tile.entrance_from_bottom or (right_tile is not None and right_tile.entrance_from_left)):
    #         raise ValueError("Entrance from left and down must feed both ways")
    #     if right_tile is None:
    #         if left_tile.feeds_right and down_tile.feeds_up:
    #             return TileType.CROSS
    #         if left_tile.feeds_right:
    #             return TileType.HORIZ
    #         return TileType.ELBOW_SE
    #     if not right_tile.entrance_from_left and not up_tile.entrance_from_bottom:
    #         raise ValueError("Middle TBD tile must feed out right or up")
    #     if (right_tile.entrance_from_left and up_tile.entrance_from_bottom) and (not left_tile.feeds_right or not down_tile.feeds_up):
    #         raise ValueError("Entrance from right and up must feed both ways")
    #     if (left_tile.feeds_right and down_tile.feeds_up) and (right_tile.entrance_from_left and up_tile.entrance_from_bottom):
    #         return TileType.CROSS
    #     if left_tile.feeds_right:
    #         if up_tile.entrance_from_bottom:
    #             return TileType.ELBOW_NW
    #         return TileType.HORIZ
    #     if up_tile.entrance_from_bottom:
    #         return TileType.VERT
    #     return TileType.ELBOW_SE

    def _get_tbd_tile(self, left_tile, right_tile, down_tile, up_tile) -> TileType:
        # bottom left, vert or elbow
        # if down_tile is None:
        #     if left_tile is None or not left_tile.feeds_right:
        #         if right_tile is None or right_tile == TileType.CROSS:
        #             if self.DEBUG:
        #                 print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")

        #             return TileType.ELBOW_SE
        #         return TileType.VERT
        #     raise ValueError("If bottom tile is not a cross, left tile cannot feed right")
        # if left_tile is None:
        #     if not down_tile.feeds_up:
        #         raise ValueError("Down tile must feed into left tile")
        #     if up_tile is None:
        #         if right_tile == TileType.EMPTY:
        #             raise ValueError("Upper left corner can only be SE elbow if not empty")
        #         if self.DEBUG:
        #             print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")

        #         return TileType.ELBOW_SE
        #     if right_tile in (TileType.EMPTY, TileType.TBD):
        #         return TileType.VERT
        #     return TileType.ELBOW_SE
        # if left_tile.feeds_right and down_tile.feeds_up:
        #     raise ValueError("Cannot have both left and down tiles feeding into TBD tile")
        # if not left_tile.feeds_right and not down_tile.feeds_up:
        #     raise ValueError("Something must feed into TBD tile")
        if up_tile is None:
            if left_tile is None or not left_tile.feeds_right:
                if self.DEBUG:
                    print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
                return TileType.ELBOW_SE
            return TileType.HORIZ
        if left_tile is None:
            if up_tile.entrance_from_bottom:
                return TileType.VERT
            if self.DEBUG:
                print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
        if left_tile.feeds_right:
            if up_tile.entrance_from_bottom:
                return TileType.ELBOW_NW
            return TileType.HORIZ
        if up_tile.entrance_from_bottom:
            return TileType.VERT
        if self.DEBUG:
            print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
        return TileType.ELBOW_SE

    DEBUG = True

    def build(self):
        """Build internal structures if needed (currently a placeholder)"""
        for col in range(self.n):
            for row in range(self.n):
                if self.grid[row, col] == TileType.TBD:
                    self.grid[row, col] = self._get_tbd_tile(
                        TileType(self.grid[row, col - 1]) if col > 0 else None,
                        TileType(self.grid[row, col + 1]) if col < self.n - 1 else None,
                        TileType(self.grid[row + 1, col]) if row < self.n - 1 else None,
                        TileType(self.grid[row - 1, col]) if row > 0 else None,
                    )
                if self.DEBUG:
                    print(f"After building cell ({row}, {col}):\n{self}")

    def __len__(self):
        """Return the size n of the n×n grid"""
        return self.n

    def __getitem__(self, key):
        """Access grid elements"""
        return self.grid[key]

    def __repr__(self):
        return f"BPD({self.grid.tolist()})"

    def __str__(self):
        """String representation of the BPD using tile symbols"""
        result = []
        for i in range(self.n):
            row = []
            for j in range(self.n):
                tile = TileType(self.grid[i, j])
                row.append(str(tile))
            result.append("".join(row))
        return "\n".join(result)

    @cached_property
    def perm(self) -> Permutation:
        """
        Compute the permutation associated with this BPD.

        The permutation is determined by following each vertical pipe from bottom to top.
        Pipes enter from the bottom (vertical) and left (horizontal).

        Returns:
            Permutation object
        """
        buildperm = []
        for col in range(self.n):
            # follow path
            current_row = self.n - 1
            current_col = col
            going_up = True
            while current_col < self.n and current_row >= 0:
                tile = TileType(self.grid[current_row, current_col])
                if tile == TileType.ELBOW_NW:
                    # Elbow NW: go left
                    assert not going_up, "Invalid pipe direction at ELBOW_NW"
                    current_row -= 1
                    going_up = True
                elif tile == TileType.ELBOW_SE:
                    # Elbow SE: go down
                    assert going_up, "Invalid pipe direction at ELBOW_SE"
                    current_col += 1
                    going_up = False
                elif not going_up:
                    current_col += 1
                else:
                    current_row -= 1
                # Empty and other tiles continue straight (down)
            buildperm.append(current_row + 1)
        return ~Permutation(buildperm)

    @property
    def permutation(self) -> Permutation:
        """Alias for perm property"""
        return self.perm

    @cached_property
    def weight(self) -> Tuple[int, ...]:
        """
        Compute the weight of this BPD.

        The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
        of empty squares (0s) in column i.

        Returns:
            Tuple of integers representing the weight
        """
        return tuple(int(np.sum(self.grid[:, j] == 0)) for j in range(self.n))

    @cached_property
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
        for row in range(self.n - 1, -1, -1):
            for col in range(self.n):
                if self.grid[row, col] == 1:
                    # Count pipes weakly northeast of this crossing
                    # Weakly northeast: r <= row and c >= col
                    pipes_northeast = self.n - col

                    # Word value is pipes_northeast - 1 (1-indexed)
                    word.append(pipes_northeast)

        return tuple(word)

    @cached_property
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
        for row in range(self.n):
            for col in range(self.n):
                if self.grid[row, col] == 1:
                    # Count pipes weakly northeast of this crossing
                    pipes_northeast = self.n - col
                    word.append(pipes_northeast)

        return tuple(word)

    def is_valid(self) -> bool:
        """
        Check if this is a valid pipe dream.

        A valid pipe dream must satisfy:
        1. Each pipe path from bottom to top does not go out of bounds
        2. The resulting permutation is valid

        Returns:
            True if valid, False otherwise
        """
        try:
            # Check that following pipes does not go out of bounds
            for start_col in range(self.n):
                current_col = start_col
                coming_from = "bottom"

                for row in range(self.n - 1, -1, -1):
                    if current_col < 0 or current_col >= self.n:
                        return False

                    tile = TileType(self.grid[row, current_col])

                    # Update position based on tile
                    if tile == TileType.CROSS and coming_from == "bottom":
                        current_col += 1
                        coming_from = "left"
                    elif tile == TileType.ELBOW_NE and coming_from == "bottom":
                        current_col += 1
                        coming_from = "left"
                    elif tile == TileType.ELBOW_NW:
                        if coming_from == "left":
                            coming_from = "bottom"
                        elif coming_from == "bottom":
                            current_col -= 1
                            coming_from = "left"
                    # Empty and other tiles continue straight

            # Check that the permutation is valid
            perm = self.perm
            return sorted(perm) == list(range(1, self.n + 1))
        except Exception:
            return False

    # @classmethod
    # def from_permutation(cls, perm: Permutation, positions: Optional[List[Tuple[int, int]]] = None):
    #     """
    #     Create a BPD from a permutation.

    #     Args:
    #         perm: Permutation object
    #         positions: Optional list of (row, col) positions where crossings should be placed.
    #                   If None, creates the essential set (minimal) BPD.

    #     Returns:
    #         BPD object
    #     """
    #     n = len(perm)
    #     grid = np.zeros((n, n), dtype=int)

    #     if positions is not None:
    #         # Place crossings at specified positions
    #         for row, col in positions:
    #             if 0 <= row < n and 0 <= col < n:
    #                 grid[row, col] = 1
    #     else:
    #         # Create essential set (Rothe diagram style)
    #         perm_list = list(perm)
    #         for i in range(n):
    #             for j in range(n):
    #                 # Place crossing if perm[i] > j+1 and i < inverse_perm[j]
    #                 inv_perm = [perm_list.index(k + 1) for k in range(n)]
    #                 if perm_list[i] > j + 1 and i < inv_perm[j]:
    #                     grid[i, j] = 1

    #     return cls(grid)

    def __eq__(self, other):
        """Check equality of two BPDs"""
        if not isinstance(other, BPD):
            return False
        return np.array_equal(self.grid, other.grid)

    def __hash__(self):
        """Hash for use in sets and dicts"""
        return hash(self.grid.tobytes())

    def copy(self):
        """Create a copy of this BPD"""
        return BPD(self.grid.copy())

    @property
    def num_crossings(self) -> int:
        """Total number of crossings in the BPD"""
        return int(np.sum(self.grid == TileType.CROSS))

    def trace_pipe(self, i, j, direction=None):
        if self.grid[i, j] == TileType.ELBOW_NW:
            if direction == "left":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self.grid[i, j] == TileType.ELBOW_SE:
            if direction == "down":
                return None
            if i == self.n - 1:
                return j + 1
            return self.trace_pipe(i + 1, j, direction="down")
        if self.grid[i, j] == TileType.HORIZ:
            if direction == "down":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self.grid[i, j] == TileType.VERT:
            if direction == "left":
                return None
            if i == self.n - 1:
                return j + 1
            return self.trace_pipe(i + 1, j, direction="down")
        if self.grid[i, j] == TileType.CROSS:
            if direction == "down":
                if i == self.n - 1:
                    return j + 1
                return self.trace_pipe(i + 1, j, direction="down")
            if direction == "left":
                return self.trace_pipe(i, j - 1, direction="left")
        raise ValueError("Invalid tile for tracing pipe")

    def undroop_elbow(self, r, c, from_spot=None):
        """
        Inverse operation of droop_elbow. Takes an NW elbow and moves it back up and left.

        Args:
            r, c: Position of the NW elbow to undroop (bottom-right corner of rectangle)
            from_spot: Optional (row, col) of where to place the SE elbow (top-left corner)

        Returns:
            New BPD with the elbow undrooped
        """

        if from_spot is None:
            # Find the closest strictly northwest EMPTY position
            # Start from (r-1, c-1) and search outward
            pot_from_spot = (r - 1, c - 1)
            while self.grid[pot_from_spot[0], pot_from_spot[1]] != TileType.EMPTY:
                if pot_from_spot[1] > 0:
                    pot_from_spot = (pot_from_spot[0], pot_from_spot[1] - 1)
                elif pot_from_spot[0] > 0:
                    pot_from_spot = (pot_from_spot[0] - 1, c - 1)
                else:
                    raise ValueError("No available empty spot northwest to undroop elbow")
            from_spot = pot_from_spot

        new_bpd = self.copy()

        # Place SE elbow at original position - removes one EMPTY
        new_bpd.grid[from_spot[0], from_spot[1]] = TileType.ELBOW_SE

        # Remove the NW elbow at destination - creates one EMPTY (net zero EMPTYs)
        new_bpd.grid[r, c] = TileType.EMPTY

        # Revert top horizontal edge from (from_spot[0], from_spot[1]+1) to (from_spot[0], c)
        for j in range(from_spot[1] + 1, c + 1):
            if new_bpd.grid[from_spot[0], j] == TileType.CROSS:
                new_bpd.grid[from_spot[0], j] = TileType.VERT
            elif new_bpd.grid[from_spot[0], j] == TileType.HORIZ:
                new_bpd.grid[from_spot[0], j] = TileType.EMPTY
            elif new_bpd.grid[from_spot[0], j] == TileType.ELBOW_SE:
                # This was a corner elbow, revert to VERT or EMPTY
                new_bpd.grid[from_spot[0], j] = TileType.EMPTY

        # Revert right vertical edge from (from_spot[0]+1, c) to (r)
        for i in range(from_spot[0] + 1, r + 1):
            if new_bpd.grid[i, c] == TileType.CROSS:
                new_bpd.grid[i, c] = TileType.HORIZ
            elif new_bpd.grid[i, c] == TileType.VERT:
                new_bpd.grid[i, c] = TileType.EMPTY

        # Revert bottom horizontal edge from (r, from_spot[1]) to (r, c)
        for j in range(from_spot[1], c):
            if new_bpd.grid[r, j] == TileType.CROSS:
                new_bpd.grid[r, j] = TileType.VERT
            elif new_bpd.grid[r, j] == TileType.HORIZ:
                new_bpd.grid[r, j] = TileType.EMPTY
            elif new_bpd.grid[r, j] == TileType.ELBOW_NW:
                # This was a corner elbow, revert to HORIZ or EMPTY
                new_bpd.grid[r, j] = TileType.EMPTY

        # Revert left vertical edge from (from_spot[0]+1, from_spot[1]) to (r-1, from_spot[1])
        for i in range(from_spot[0] + 1, r):
            if new_bpd.grid[i, from_spot[1]] == TileType.CROSS:
                new_bpd.grid[i, from_spot[1]] = TileType.HORIZ
            elif new_bpd.grid[i, from_spot[1]] == TileType.VERT:
                new_bpd.grid[i, from_spot[1]] = TileType.EMPTY

        return new_bpd

    def droop_elbow(self, r, c, to_spot=None):
        if self.grid[r, c] != TileType.ELBOW_SE:
            raise ValueError("No SE elbow at given position to droop")
        if to_spot is None:
            # find lowest empty spot southeast of (r, c)
            pot_to_spot = (r + 1, c + 1)
            while self.grid[pot_to_spot[0], pot_to_spot[1]] != TileType.EMPTY:
                if pot_to_spot[1] < self.n - 1:
                    pot_to_spot = (pot_to_spot[0], pot_to_spot[1] + 1)
                elif pot_to_spot[0] < self.n - 1:
                    pot_to_spot = (pot_to_spot[0] + 1, c + 1)
                else:
                    raise ValueError("No available spot to droop elbow")
            to_spot = pot_to_spot

        new_bpd = self.copy()

        # Remove SE elbow at origin - creates one EMPTY
        new_bpd.grid[r, c] = TileType.EMPTY

        # Place destination elbow - removes one EMPTY (net zero EMPTYs)
        new_bpd.grid[to_spot[0], to_spot[1]] = TileType.ELBOW_NW

        # Top horizontal edge from (r, c+1) to (r, to_spot[1])
        for j in range(c + 1, to_spot[1] + 1):
            if new_bpd.grid[r, j] == TileType.VERT:
                new_bpd.grid[r, j] = TileType.CROSS
            elif new_bpd.grid[r, j] == TileType.EMPTY:
                new_bpd.grid[r, j] = TileType.HORIZ

        # Corner: turn the horizontal pipe down at (r, to_spot[1])
        if new_bpd.grid[r, to_spot[1]] == TileType.HORIZ:
            new_bpd.grid[r, to_spot[1]] = TileType.ELBOW_SE
        elif new_bpd.grid[r, to_spot[1]] == TileType.CROSS:
            # If there was already a vertical pipe, keep it as CROSS
            pass

        # Right vertical edge from (r+1, to_spot[1]) to (to_spot[0])
        for i in range(r + 1, to_spot[0] + 1):
            if new_bpd.grid[i, to_spot[1]] == TileType.HORIZ:
                new_bpd.grid[i, to_spot[1]] = TileType.CROSS
            elif new_bpd.grid[i, to_spot[1]] == TileType.EMPTY:
                new_bpd.grid[i, to_spot[1]] = TileType.VERT

        # Bottom horizontal edge from (to_spot[0], c) to (to_spot[0], to_spot[1])
        for j in range(c, to_spot[1] + 1):
            if new_bpd.grid[to_spot[0], j] == TileType.VERT:
                new_bpd.grid[to_spot[0], j] = TileType.CROSS
            elif new_bpd.grid[to_spot[0], j] == TileType.EMPTY:
                new_bpd.grid[to_spot[0], j] = TileType.HORIZ

        # Corner: turn the vertical pipe left at (to_spot[0], c)
        if new_bpd.grid[to_spot[0], c] == TileType.VERT:
            new_bpd.grid[to_spot[0], c] = TileType.ELBOW_NW
        elif new_bpd.grid[to_spot[0], c] == TileType.CROSS:
            # If there was already a horizontal pipe, keep it as CROSS
            pass

        # Left vertical edge from (r+1, c) to (to_spot[0])
        for i in range(r + 1, to_spot[0]):
            if new_bpd.grid[i, c] == TileType.HORIZ:
                new_bpd.grid[i, c] = TileType.CROSS
            elif new_bpd.grid[i, c] == TileType.EMPTY:
                new_bpd.grid[i, c] = TileType.VERT

        return new_bpd

        return new_bpd

    # def delta_op(self):
    #     r = min(i for i in range(self.n) if any(self.grid[i, j] == TileType.EMPTY for j in range(self.n)))
    #     mark = max(j for j in range(self.n) if self.grid[r, j] == TileType.EMPTY)
    #     new_bpd = self.copy()
    #     while True:
    #         assert new_bpd.grid[r, mark] == TileType.EMPTY
    #         while mark < self.n - 1 and new_bpd.grid[r, mark + 1] == TileType.EMPTY:
    #             mark = mark + 1
    #         pipe = new_bpd.trace_pipe(r, mark + 1)
    #         print(f"da pipe iz {pipe=}")
    #         if pipe != mark + 2:
    #             # find the elbow
    #             print("I iz to elbow the the")
    #             elbow_pos = r + 1
    #             while elbow_pos < self.n and new_bpd.grid[elbow_pos, mark + 1] != TileType.ELBOW_NW:
    #                 elbow_pos += 1
    #             if elbow_pos == self.n:
    #                 raise ValueError("No elbow found for delta operation")
    #             for middle_spot in range(r + 1, elbow_pos):
    #                 crossing_pipe = new_bpd.trace_pipe(middle_spot, mark + 1, direction="left")
    #                 if new_bpd.grid[middle_spot, mark] == TileType.ELBOW_SE:
    #                     # want to swap crossing_pipe
    #                     for other_elbow_pos in range(middle_spot + 1, self.n):
    #                         if new_bpd.grid[other_elbow_pos, mark] == TileType.ELBOW_NW:
    #                             # found it
    #                             new_bpd = self.copy()
    #                             new_bpd = new_bpd.droop_elbow(middle_spot, mark, to_spot=(other_elbow_pos, mark+1))
    #                             new_bpd.grid[middle_spot, mark + 1] = TileType.RC_EMPTY
    #             new_bpd = new_bpd.undroop_elbow(elbow_pos, mark + 1, from_spot=(r, mark))
    #             mark = mark + 1
    #             r = elbow_pos

    #         else:
    #             print("I ain't sticknk")
    #             cross_pos = r + 1
    #             while cross_pos < self.n and new_bpd.grid[cross_pos, mark + 1] != TileType.CROSS:
    #                 cross_pos += 1
    #             if cross_pos == self.n:
    #                 raise ValueError(f"No crossing found for delta operation \n{new_bpd}\n{r=} {mark=} ")
    #             new_bpd = new_bpd.undroop_elbow(cross_pos, mark + 1, from_spot=(r, mark))
    #             new_bpd.grid[cross_pos, mark + 1] = TileType.ELBOW_SE
    #             return new_bpd

    def to_rc_graph(self):
        """
        Convert this BPD to an RC-graph representation.

        Returns:
            RCGraph object (if available in the module)
        """
        from schubmult.schub_lib.rc_graph import RCGraph

        # RC-graph rows contain the column positions of crossings in each row
        rows = []
        for i in range(self.n):
            row_crossings = []
            for j in range(self.n):
                if self.grid[i, j] == TileType.CROSS:
                    row_crossings.append(j + 1)  # 1-indexed
            rows.append(tuple(row_crossings))

        return RCGraph(tuple(rows))

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

    def __str__(self):
        symbols = {TileType.EMPTY: "·", TileType.CROSS: "┼", TileType.ELBOW_NW: "╯", TileType.ELBOW_SE: "╭", TileType.HORIZ: "─", TileType.VERT: "│", TileType.TBD: "?"}
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
            return TileType.ELBOW_SE
        if left_tile.feeds_right:
            if up_tile.entrance_from_bottom:
                return TileType.ELBOW_NW
            return TileType.HORIZ
        if up_tile.entrance_from_bottom:
            return TileType.VERT
        if self.DEBUG:
            print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
        return TileType.ELBOW_SE

    DEBUG = False

    def build(self):
        """Build internal structures if needed (currently a placeholder)"""
        for col in range(self.n):
            for row in range(self.n):
                if self[row, col] == TileType.TBD:
                    self.grid[row, col] = self._get_tbd_tile(
                        TileType(self[row, col - 1]) if col > 0 else None,
                        TileType(self[row, col + 1]) if col < self.n - 1 else None,
                        TileType(self[row + 1, col]) if row < self.n - 1 else None,
                        TileType(self[row - 1, col]) if row > 0 else None,
                    )
                if self.DEBUG:
                    print(f"After building cell ({row}, {col}):\n{self}")
        assert self.is_valid(), f"Built BPD is not valid: \n{self}"

    def __len__(self):
        """Return the size n of the n×n grid"""
        return self.n

    def __getitem__(self, key):
        """Access grid elements, casting to TileType"""
        result = self.grid[key]
        # If it's a numpy array (from slicing), cast each element
        if isinstance(result, np.ndarray):
            return result.astype(TileType)
        # Otherwise it's a scalar, cast directly
        return TileType(result)

    def __repr__(self):
        return f"BPD({self.grid.tolist()})"

    def __str__(self):
        """String representation of the BPD using tile symbols"""
        result = []
        for i in range(self.n):
            row = []
            for j in range(self.n):
                tile = TileType(self[i, j])
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
                tile = TileType(self[current_row, current_col])
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
        return tuple(int(np.sum(self[:, j] == 0)) for j in range(self.n))

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
                if self[row, col] == 1:
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
                if self[row, col] == 1:
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
        # sanity checks
        for row in range(self.n):
            if self[row, 0] in (TileType.ELBOW_NW, TileType.CROSS, TileType.HORIZ):
                return False
            if self[row, self.n - 1] in (TileType.ELBOW_NW, TileType.VERT):
                return False
        for col in range(self.n):
            if self[self.n - 1, col] in (TileType.ELBOW_NW, TileType.HORIZ):
                return False
            if self[0, col] in (TileType.ELBOW_NW, TileType.VERT, TileType.CROSS):
                return False
        for row in range(1, self.n - 1):
            for col in range(1, self.n - 1):
                if self[row, col].feeds_right and not self[row, col + 1].entrance_from_left:
                    return False
                if self[row, col].feeds_up and not self[row - 1, col].entrance_from_bottom:
                    return False
                if self[row, col].entrance_from_left and not self[row, col - 1].feeds_right:
                    return False
                if self[row, col].entrance_from_bottom and not self[row + 1, col].feeds_up:
                    return False
        return True
        # try:
        #     # Check that following pipes does not go out of bounds
        #     for start_col in range(self.n):
        #         current_col = start_col
        #         coming_from = "bottom"

        #         for row in range(self.n - 1, -1, -1):
        #             if current_col < 0 or current_col >= self.n:
        #                 return False

        #             tile = TileType(self[row, current_col])

        #             # Update position based on tile
        #             if tile == TileType.CROSS and coming_from == "bottom":
        #                 current_col += 1
        #                 coming_from = "left"
        #             elif tile == TileType.ELBOW_NE and coming_from == "bottom":
        #                 current_col += 1
        #                 coming_from = "left"
        #             elif tile == TileType.ELBOW_NW:
        #                 if coming_from == "left":
        #                     coming_from = "bottom"
        #                 elif coming_from == "bottom":
        #                     current_col -= 1
        #                     coming_from = "left"
        #             # Empty and other tiles continue straight

        #     # Check that the permutation is valid
        #     perm = self.perm
        #     return sorted(perm) == list(range(1, self.n + 1))
        # except Exception:
        #     return False

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
        if self[i, j] == TileType.ELBOW_NW:
            if direction == "left":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self[i, j] == TileType.ELBOW_SE:
            if direction == "down":
                return None
            if i == self.n - 1:
                return j + 1
            return self.trace_pipe(i + 1, j, direction="down")
        if self[i, j] == TileType.HORIZ:
            if direction == "down":
                return None
            return self.trace_pipe(i, j - 1, direction="left")
        if self[i, j] == TileType.VERT:
            if direction == "left":
                return None
            if i == self.n - 1:
                return j + 1
            return self.trace_pipe(i + 1, j, direction="down")
        if self[i, j] == TileType.CROSS:
            if direction == "down":
                if i == self.n - 1:
                    return j + 1
                return self.trace_pipe(i + 1, j, direction="down")
            if direction == "left":
                return self.trace_pipe(i, j - 1, direction="left")
            raise ValueError("Must specify direction when tracing through a crossing")
        raise ValueError(f"Invalid tile for tracing pipe at ({i}, {j}): {self[i, j]}")

    def droop_elbow_to_empty(self, r, c, to_spot=None):
        # assume (r, c) is empty
        # to_spot is empty
        if to_spot is None:
            # find cloest se elbow spot nw of (r, c)
            pot_to_spot = (r - 1, c - 1)
            while self[pot_to_spot[0], pot_to_spot[1]] != TileType.ELBOW_SE:
                if pot_to_spot[1] > 0:
                    pot_to_spot = (pot_to_spot[0], pot_to_spot[1] - 1)
                elif pot_to_spot[0] > 0:
                    pot_to_spot = (pot_to_spot[0] - 1, c - 1)
                else:
                    raise ValueError("No available spot to droop elbow")
            to_spot = pot_to_spot

        new_bpd = np.full(self.grid.shape, TileType.TBD, dtype=TileType)

        # move the empty

        for i in range(self.n):
            for j in range(self.n):
                if i <= r and (j == to_spot[1]) and i >= to_spot[0]:
                    continue
                if j <= c and (i == to_spot[0]) and j >= to_spot[1]:
                    continue
                new_bpd[i, j] = self[i, j] if self[i, j] in (TileType.CROSS, TileType.EMPTY) else TileType.TBD

        for j in range(to_spot[1] + 1, c):
            if self[r, j] == TileType.VERT:
                new_bpd[r, j] = TileType.CROSS
        for i in range(to_spot[0] + 1, r):
            if self[i, c] == TileType.HORIZ:
                new_bpd[i, c] = TileType.CROSS

        new_bpd[to_spot[0], to_spot[1]] = TileType.EMPTY

        new_bpd[r, c] = TileType.TBD
        return BPD(new_bpd)

    def undroop_elbow_to_empty(self, r, c, from_spot=None):
        """
        Reverse of droop_elbow_to_empty.
        Assumes (r, c) is empty and finds SE elbow strictly southeast to move to (r, c).
        Crossings may appear on west and north edges (opposite of droop).
        """
        if from_spot is None:
            # find closest SE elbow strictly southeast of (r, c)
            pot_from_spot = (r + 1, c + 1)
            while self[pot_from_spot[0], pot_from_spot[1]] != TileType.ELBOW_NW:
                if pot_from_spot[1] < self.n - 1:
                    pot_from_spot = (pot_from_spot[0], pot_from_spot[1] + 1)
                elif pot_from_spot[0] < self.n - 1:
                    pot_from_spot = (pot_from_spot[0] + 1, c + 1)
                else:
                    raise ValueError("No available NW elbow southeast to undroop")
            from_spot = pot_from_spot

        new_bpd = np.full(self.grid.shape, TileType.TBD, dtype=TileType)

        # Copy CROSS and EMPTY tiles, excluding the rectangle edges
        for i in range(self.n):
            for j in range(self.n):
                if i >= r and (j == from_spot[1]) and i <= from_spot[0]:
                    continue
                if j >= c and (i == from_spot[0]) and j <= from_spot[1]:
                    continue
                new_bpd[i, j] = self[i, j] if self[i, j] in (TileType.CROSS, TileType.EMPTY) else TileType.TBD

        # Add crossings on west edge where vertical pipes cross new horizontal path
        for i in range(r + 1, from_spot[0]):
            if self[i, c] == TileType.HORIZ:
                new_bpd[i, c] = TileType.CROSS

        # Add crossings on north edge where horizontal pipes cross new vertical path
        for j in range(c + 1, from_spot[1]):
            if self[r, j] == TileType.VERT:
                new_bpd[r, j] = TileType.CROSS

        # Move the empty from from_spot to (r, c)
        new_bpd[from_spot[0], from_spot[1]] = TileType.EMPTY
        new_bpd[r, c] = TileType.TBD

        return BPD(new_bpd)

    def delta_op(self):
        """
        Perform the delta operation on a pipe dream according to Definition 3.1.
        Only manipulates CROSS and EMPTY tiles, uses build() to reconstruct pipes.

        Returns:
            Tuple of (new_bpd, position) where position is 1-indexed (col, row)
        """
        # Find the smallest row index r that contains empty tiles
        r = min(i for i in range(self.n) if any(self[i, j] == TileType.EMPTY for j in range(self.n)))

        # Find the rightmost empty in that row - this is our initial mark
        y = max(j for j in range(self.n) if self[r, j] == TileType.EMPTY)
        x = r
        orig_r = r

        if self.DEBUG:
            print(f"Initial mark at ({x}, {y}) (0-indexed)")
            print(f"Initial state:\n{self}\n")

        # Work on a copy
        current_bpd = self.copy()

        while True:
            # Step (1): Move mark to rightmost empty in contiguous block
            while y + 1 < self.n and current_bpd.grid[x, y + 1] == TileType.EMPTY:
                y = y + 1
                if self.DEBUG:
                    print(f"Step (1): Moving mark right to ({x}, {y})")

            if self.DEBUG:
                print(f"Mark is now at ({x}, {y})")
                print(f"Current state:\n{current_bpd}\n")

            p = current_bpd.trace_pipe(x, y + 1, direction="down")
            if p is None:
                p = current_bpd.trace_pipe(x, y + 1, direction="left")
            if self.DEBUG:
                print(f"Traced pipe from ({x}, {y + 1}) to column {p} (1-indexed)")
            if p != y + 2:
                # find x_prime
                x_prime = x
                while x_prime < self.n and current_bpd.grid[x_prime, y + 1] != TileType.ELBOW_NW:
                    x_prime += 1

                if self.DEBUG:
                    print(f"Step (2a): Found x' at ({x_prime}, {y + 1})")
                new_bpd = current_bpd.copy()
                for z in range(x + 1, x_prime):
                    if new_bpd.grid[z, y + 1] == TileType.CROSS:
                        new_bpd.grid[z, y] = TileType.CROSS
                        new_bpd.grid[z, y + 1] = TileType.TBD
                new_bpd.grid[x, y] = TileType.TBD
                new_bpd.grid[x_prime, y + 1] = TileType.EMPTY
                new_bpd.rebuild()
                current_bpd = new_bpd
                x = x_prime
                y = y + 1
            else:
                if self.DEBUG:
                    print(f"Step (2b): Performing final move at ({x}, {y})")
                # x, y may not be the empty spot
                empty_row = x
                empty_col = y

                current_bpd.grid[empty_row, empty_col] = TileType.TBD
                x_prime = x + 1
                for z in range(x + 1, self.n):
                    if current_bpd.grid[z, y + 1] == TileType.CROSS:
                        x_prime = z
                        break
                if self.DEBUG:
                    print(f"Final x' found at ({x_prime}, {y + 1})")
                current_bpd.grid[x_prime, y + 1] = TileType.TBD
                for z in range(x + 1, x_prime):
                    if current_bpd.grid[z, y + 1] == TileType.CROSS:
                        current_bpd.grid[z, y] = TileType.CROSS
                        current_bpd.grid[z, y + 1] = TileType.TBD
                current_bpd.rebuild()

                return current_bpd, (y + 1, orig_r + 1)  # return 1-indexed position

    def rebuild(self):
        """Rebuild the BPD to resolve any TBD tiles"""
        for i in range(self.n):
            for j in range(self.n):
                if self[i, j] not in (TileType.EMPTY, TileType.CROSS):
                    self.grid[i, j] = TileType.TBD
        self.build()

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
                if self[i, j] == TileType.CROSS:
                    row_crossings.append(j + 1)  # 1-indexed
            rows.append(tuple(row_crossings))

        return RCGraph(tuple(rows))

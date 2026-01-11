"""
Bumpless Pipe Dreams (BPD) module

A bumpless pipe dream is a pipe dream diagram with no elbows - only crossings and empty boxes.
Represented as an n x n grid where:
- 1 = crossing (pipes cross)
- 0 = empty box (pipes go straight through)

For general pipe dreams, there are 6 possible tile types.
"""

from enum import IntEnum
from typing import Tuple

import numpy as np

from schubmult.schub_lib.perm_lib import Permutation


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

    def __str__(self):
        symbols = {TileType.BLANK: "·", TileType.CROSS: "┼", TileType.ELBOW_NW: "╯", TileType.ELBOW_SE: "╭", TileType.HORIZ: "─", TileType.VERT: "│", TileType.TBD: "?"}
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


class BPD:
    """
    Bumpless Pipe Dream representation.

    A bumpless pipe dream is an n×n grid where:
    - TileType.CROSS (1) represents a crossing
    - TileType.BLANK (0) represents an empty box (pipes go straight)
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


    @classmethod
    def _get_tbd_tile(cls, left_tile, up_tile) -> TileType:
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
                return TileType.ELBOW_NW
            return TileType.HORIZ
        if up_tile.entrance_from_bottom:
            return TileType.VERT
        # if self.DEBUG:
        #     print(f"Returning SE elbow for cell with {left_tile=}, {down_tile=}, {up_tile=}, {right_tile=}, line={__import__('inspect').currentframe().f_back.f_lineno}")
        return TileType.ELBOW_SE

    DEBUG = False

    def build(self, validate=False):
        """Build internal structures if needed (currently a placeholder)"""
        for col in range(self.n):
            for row in range(self.n):
                if self[row, col] == TileType.TBD:
                    self.grid[row, col] = self._get_tbd_tile(
                        TileType(self[row, col - 1]) if col > 0 else None,
                        TileType(self[row - 1, col]) if row > 0 else None,
                    )
                # if self.DEBUG:
                #     print(f"After building cell ({row}, {col}):\n{self}")
        if validate:
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

    @classmethod
    def from_asm(cls, asm):
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
    def rothe_bpd(cls, perm):
        n = len(perm)
        grid = np.full((n, n), fill_value=TileType.TBD, dtype=TileType)
        bpd = BPD(grid)
        for a, b in perm.diagram:
            bpd.grid[a - 1, b - 1] = TileType.BLANK
        graph = perm.graph
        for i in range(n):
            for j in range(n):
                if bpd[i , j] != TileType.BLANK:
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
        return tuple(int(np.sum(self[:, j] == TileType.BLANK)) for j in range(self.n))

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
        for row in range(self.n - 1, -1, -1):
            for col in range(self.n):
                if self[row, col] == 1:
                    # Count pipes weakly northeast of this crossing
                    # Weakly northeast: r <= row and c >= col
                    pipes_northeast = self.n - col

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

    def droop_elbow_to_empty(self, r, c, to_spot=None, override_tile=TileType.BLANK):
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

        new_bpd = self.copy()

        # move the empty

        for i in range(self.n):
            for j in range(self.n):
                if i <= r and (j == to_spot[1]) and i >= to_spot[0]:
                    new_bpd.grid[i, j] = TileType.TBD
                if j <= c and (i == to_spot[0]) and j >= to_spot[1]:
                    new_bpd.grid[i, j] = TileType.TBD
                #new_bpd[i, j] = self[i, j] if self[i, j] in (TileType.CROSS, TileType.BLANK) else TileType.TBD

        for j in range(to_spot[1] + 1, c):
            if self[r, j] == TileType.VERT:
                new_bpd[r, j] = TileType.CROSS
        for i in range(to_spot[0] + 1, r):
            if self[i, c] == TileType.HORIZ:
                new_bpd[i, c] = TileType.CROSS

        new_bpd.grid[to_spot[0], to_spot[1]] = override_tile

        new_bpd.grid[r, c] = TileType.TBD
        new_bpd.rebuild()
        return new_bpd

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

        new_bpd = self.copy()

        # Copy CROSS and EMPTY tiles, excluding the rectangle edges
        for i in range(self.n):
            for j in range(self.n):
                if i >= r and (j == from_spot[1]) and i <= from_spot[0]:
                    new_bpd.grid[i, j] = TileType.TBD
                if j >= c and (i == from_spot[0]) and j <= from_spot[1]:
                    new_bpd.grid[i, j] = TileType.TBD
                #new_bpd[i, j] = self[i, j] if self[i, j] in (TileType.CROSS, TileType.BLANK) else TileType.TBD

        # Add crossings on west edge where vertical pipes cross new horizontal path
        for i in range(r + 1, from_spot[0]):
            if self[i, c] == TileType.HORIZ:
                new_bpd.grid[i, c] = TileType.CROSS

        # Add crossings on north edge where horizontal pipes cross new vertical path
        for j in range(c + 1, from_spot[1]):
            if self[r, j] == TileType.VERT:
                new_bpd.grid[r, j] = TileType.CROSS

        # Move the empty from from_spot to (r, c)
        new_bpd.grid[from_spot[0], from_spot[1]] = TileType.BLANK
        new_bpd.grid[r, c] = TileType.TBD
        new_bpd.rebuild()
        return new_bpd

    def delta_op(self):
        """
        Perform the delta operation on a pipe dream according to Definition 3.1.
        Only manipulates CROSS and EMPTY tiles, uses build() to reconstruct pipes.

        Returns:
            Tuple of (new_bpd, position) where position is 1-indexed (col, row)
        """
        # Find the smallest row index r that contains empty tiles
        r = min(i for i in range(self.n) if any(self[i, j] == TileType.BLANK for j in range(self.n)))

        # Find the rightmost empty in that row - this is our initial mark
        y = max(j for j in range(self.n) if self[r, j] == TileType.BLANK)
        x = r
        orig_r = r

        if self.DEBUG:
            print(f"Initial mark at ({x}, {y}) (0-indexed)")
            print(f"Initial state:\n{self}\n")

        # Work on a copy
        current_bpd = self.copy()
        final_col = None  # Will track the column for return value

        while True:
            # Step (1): Move mark to rightmost empty in contiguous block
            if self.DEBUG:
                print("\n=== Step (1): Moving mark to rightmost empty in contiguous block ===")
                print(f"Starting position: ({x}, {y})")
                print(f"Current state:\n{current_bpd}\n")

            while y + 1 < self.n and current_bpd.grid[x, y + 1] == TileType.BLANK:
                y = y + 1
                if self.DEBUG:
                    print(f"Step (1): Moving mark right to ({x}, {y})")
                    print(f"Current state:\n{current_bpd}\n")
            assert current_bpd.grid[x, y] == TileType.BLANK, "Mark must be on an empty tile after Step (1)"
            if self.DEBUG:
                print(f"Mark is now at ({x}, {y})")
                print(f"Current state:\n{current_bpd}\n")

            # Save this y value - it's the left column of a potential rectangle move
            if final_col is None:
                final_col = y

            if self.DEBUG:
                print(f"Tracing pipe from ({x}, {y + 1})...")
            p = current_bpd.trace_pipe(x, y + 1, direction="down")
            if p is None:
                p = current_bpd.trace_pipe(x, y + 1, direction="left")
            if self.DEBUG:
                print(f"Traced pipe from ({x}, {y + 1}) to column {p} (1-indexed)")
                print(f"Checking if p ({p}) != y + 2 ({y + 2})")
            if p != y + 2:
                # Step 2: Column move
                if self.DEBUG:
                    print("\n=== Step (2): Column move (p != y+2) ===")
                    print(f"Pipe p = {p}, looking for where it has EMPTY in column {y + 1}")

                # Find x' where pipe p has its EMPTY (⊡-tile) in column y+1
                # This is where the elbow is after rebuild
                temp_bpd = current_bpd.copy()
                temp_bpd.rebuild()
                x_prime = x
                while x_prime < self.n and temp_bpd[x_prime, y + 1] != TileType.ELBOW_NW:
                    x_prime += 1

                if self.DEBUG:
                    print(f"Found x' = {x_prime} (pipe p has elbow at ({x_prime}, {y + 1}))")
                    print(f"Column move rectangle: NW=({x}, {y}), SE=({x_prime}, {y + 1})")
                    print(f"Before column move:\n{current_bpd}\n")

                new_bpd = current_bpd.copy()

                # Step 2a: Move tiles in the column move rectangle
                # - CROSSes in column y+1 move left and down
                # - HORIZ in column y become CROSS
                # - BLANKs in column y move right and down
                if self.DEBUG:
                    print("Step 2a: Moving tiles in column move rectangle")
                    print(f"  Rectangle: rows ({x}, {x_prime}), columns ({y}, {y + 1})")

                # First, move CROSSes from column y+1 to column y (left and down one row)
                for z in range(x + 1, x_prime):
                    if current_bpd[z, y + 1] == TileType.CROSS:
                        if self.DEBUG:
                            print(f"  CROSS: ({z}, {y + 1}) → ({z + 1}, {y})")
                        new_bpd.grid[z, y + 1] = TileType.TBD
                        new_bpd.grid[z + 1, y] = TileType.CROSS

                # Second, HORIZ tiles in column y become CROSS (pipe from left now crosses vertical pipe)
                for z in range(x, x_prime + 1):
                    if current_bpd[z, y] == TileType.HORIZ:
                        if self.DEBUG:
                            print(f"  HORIZ → CROSS: ({z}, {y})")
                        new_bpd.grid[z, y] = TileType.CROSS

                # Third, move BLANKs from column y to column y+1
                # They need to move to maintain pipe structure
                # But exclude the marked empty at (x, y) - that's handled in step 2b
                for z in range(x + 1, x_prime):  # Start from x+1, not x
                    if current_bpd[z, y] == TileType.BLANK:
                        if self.DEBUG:
                            print(f"  BLANK: ({z}, {y}) → ({z + 1}, {y + 1})")
                        new_bpd.grid[z, y] = TileType.TBD
                        new_bpd.grid[z + 1, y + 1] = TileType.BLANK

                # Step 2b: Undroop pipe p from (x', y+1) into (x, y)
                # Move the marked empty from (x, y) to (x', y+1)
                if self.DEBUG:
                    print(f"Step 2b: Undroop pipe p from ({x_prime}, {y + 1}) to ({x}, {y})")

                new_bpd.grid[x, y] = TileType.TBD
                new_bpd.grid[x_prime, y + 1] = TileType.BLANK

                if self.DEBUG:
                    print(f"Before rebuild:\n{new_bpd}\n")
                    print("Calling rebuild()...")

                new_bpd.rebuild()

                if self.DEBUG:
                    print(f"After rebuild:\n{new_bpd}\n")
                    print("Checking validity after step 2...")

                if not new_bpd.is_valid():
                    raise ValueError(f"BPD is invalid after step 2 rebuild:\n{new_bpd}")

                if self.DEBUG:
                    print(f"Valid! Updating position: x = {x_prime}, y = {y + 1}")

                current_bpd = new_bpd
                x = x_prime
                y = y + 1
            else:
                # Step 3: Final move (p == y+2)
                if self.DEBUG:
                    print("\n=== Step (3): Final move (p == y+2) ===")
                    print(f"Current position: ({x}, {y}), orig_r = {orig_r}")
                    print(f"Will return position: ({y + 1}, {orig_r + 1}) (1-indexed)")
                    print(f"Current state:\n{current_bpd}\n")

                # Find where pipes y and y+1 intersect in column y+1
                if self.DEBUG:
                    print(f"Looking for CROSS in column {y + 1} (where pipes {y} and {y+1} intersect)")
                    print(f"Crossings in column {y + 1} BEFORE any shifts:")
                    for z in range(self.n):
                        if current_bpd[z, y + 1] == TileType.CROSS:
                            print(f"  Row {z}: CROSS")

                x_prime = None
                for z in range(x + 1, self.n):
                    if current_bpd[z, y + 1] == TileType.CROSS:
                        x_prime = z
                        break

                if self.DEBUG:
                    print(f"Found crossing at ({x_prime}, {y + 1})")
                    print(f"Shifting kinks right: CROSSes and BLANKs from column {y} to column {y + 1} in rows ({x + 1}, {x_prime})")

                new_bpd = current_bpd.copy()

                # Shift crossings right: column y → column y+1 for rows between x and x'
                for z in range(x + 1, x_prime):
                    if current_bpd[z, y] == TileType.CROSS:
                        if self.DEBUG:
                            print(f"  Moving CROSS from ({z}, {y}) to ({z}, {y + 1})")
                        new_bpd.grid[z, y] = TileType.TBD
                        new_bpd.grid[z, y + 1] = TileType.CROSS
                    elif current_bpd[z, y] == TileType.BLANK:
                        if self.DEBUG:
                            print(f"  Moving BLANK from ({z}, {y}) to ({z}, {y + 1})")
                        new_bpd.grid[z, y] = TileType.TBD
                        new_bpd.grid[z, y + 1] = TileType.BLANK

                # Remove the crossing at (x', y+1) and the empty at (x, y)
                if self.DEBUG:
                    print(f"Removing CROSS at ({x_prime}, {y + 1}) and EMPTY at ({x}, {y})")

                new_bpd.grid[x_prime, y + 1] = TileType.TBD
                new_bpd.grid[x, y] = TileType.TBD

                # Per Figure 3: shift ALL remaining crossings from column y+1 to column y
                # Check against current_bpd (before shifts) to see what's in column y+1
                if self.DEBUG:
                    print(f"Shifting all CROSSes from column {y + 1} to column {y} (per Figure 3)")
                for z in range(self.n):
                    # Skip the crossing we just removed
                    if z == x_prime:
                        continue
                    # Look at what was in current_bpd before any shifts
                    if current_bpd[z, y + 1] == TileType.CROSS:
                        # Check if this crossing was already shifted right in the kink shift
                        # If so, it's now at a different location in new_bpd
                        if x + 1 <= z < x_prime:
                            # This one was already shifted from column y to y+1, so it's still there
                            # We need to move it back to column y
                            if self.DEBUG:
                                print(f"  Moving CROSS from ({z}, {y + 1}) back to ({z}, {y})")
                            new_bpd.grid[z, y + 1] = TileType.TBD
                            new_bpd.grid[z, y] = TileType.CROSS
                        else:
                            # This crossing was originally in column y+1 and not touched by kink shift
                            if self.DEBUG:
                                print(f"  Moving CROSS from ({z}, {y + 1}) to ({z}, {y})")
                            # Special case: for bottom row, mark both as TBD and let rebuild decide
                            if z == self.n - 1:
                                new_bpd.grid[z, y + 1] = TileType.TBD
                                new_bpd.grid[z - 1, y] = TileType.CROSS
                            else:
                                new_bpd.grid[z, y + 1] = TileType.TBD
                                new_bpd.grid[z - 1, y] = TileType.CROSS

                if self.DEBUG:
                    print(f"Before rebuild:\n{new_bpd}\n")
                    print("Tile inventory before rebuild():")
                    print(f"  Column {y}:")
                    for z in range(self.n):
                        tile = new_bpd.grid[z, y]
                        if tile == TileType.CROSS:
                            print(f"    Row {z}: CROSS")
                        elif tile == TileType.BLANK:
                            print(f"    Row {z}: BLANK")
                        elif tile == TileType.TBD:
                            print(f"    Row {z}: TBD")
                    print(f"  Column {y + 1}:")
                    for z in range(self.n):
                        tile = new_bpd.grid[z, y + 1]
                        if tile == TileType.CROSS:
                            print(f"    Row {z}: CROSS")
                        elif tile == TileType.BLANK:
                            print(f"    Row {z}: BLANK")
                        elif tile == TileType.TBD:
                            print(f"    Row {z}: TBD")
                    print("Calling rebuild()...")

                new_bpd.rebuild()

                if self.DEBUG:
                    print(f"After rebuild:\n{new_bpd}\n")
                    print("Checking validity after step 3...")

                if not new_bpd.is_valid():
                    raise ValueError(f"BPD is invalid after step 3 rebuild:\n{new_bpd}")

                # Check that the resulting permutation is s_y * original_perm
                from schubmult.schub_lib.perm_lib import Permutation
                s_y = Permutation.ref_product(y + 1)
                expected_perm = s_y * self.perm
                if new_bpd.perm != expected_perm:
                    if self.DEBUG:
                        print("WARNING: Permutation mismatch!")
                        print(f"  Expected: {expected_perm} (s_{y+1} * {self.perm})")
                        print(f"  Got: {new_bpd.perm}")

                # The return column should be where pipe orig_r exits in the new BPD
                # This is given by new_bpd.perm[orig_r] - 1 (converting to 0-indexed)
                exit_col = new_bpd.perm[orig_r] - 1

                if self.DEBUG:
                    print(f"Pipe {orig_r} (1-indexed: {orig_r + 1}) exits at column {exit_col} (1-indexed: {exit_col + 1})")
                    print(f"Valid! Returning position: ({exit_col + 1}, {orig_r + 1}) (1-indexed)")

                return new_bpd, (y + 1, orig_r + 1)

    def rebuild(self):
        """Rebuild the BPD to resolve any TBD tiles"""
        for i in range(self.n):
            for j in range(self.n):
                if self[i, j] not in (TileType.BLANK, TileType.CROSS):
                    self.grid[i, j] = TileType.TBD
        self.build()

    def zero_out_last_row(self):
        a, b = self.perm.maximal_corner
        assert self[a - 1, b - 1] == TileType.ELBOW_NW
        # find nearest CROSS strictly SE of this
        r, c = a, b
        while r < self.n and c < self.n and self[r, c] != TileType.CROSS:
            if c < self.n - 1:
                c += 1
            else:
                r += 1
                c = b + 1
        if r == self.n or c == self.n:
            raise ValueError("No CROSS found strictly SE of maximal corner")
        res_bpd = self.copy()
        res_bpd.grid[a - 1, b - 1] = TileType.CROSS
        res_bpd.grid[r, c] = TileType.TBD
        res_bpd.rebuild()
        return res_bpd

    def to_rc_graph(self):
        """
        Convert this BPD to an RC-graph representation.

        Returns:
            RCGraph object (if available in the module)
        """
        from schubmult.schub_lib.rc_graph import RCGraph

        # RC-graph rows contain the column positions of crossings in each row
        code_len = len(self.perm.trimcode)
        rows = [[]] * code_len
        work_bpd = self

        while work_bpd.perm.inv > 0:
            work_bpd, (reflection, row) = work_bpd.delta_op()
            rows[row - 1] = rows[row - 1] + [reflection]

        return RCGraph([tuple(r) for r in rows])

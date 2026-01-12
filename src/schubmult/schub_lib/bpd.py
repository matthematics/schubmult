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
    def all_bpds(cls, w):
        pipes = set()
        new_pipes = [BPD.rothe_bpd(w)]

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

    def all_se_elbow_spots(self):
        elbows = set()
        for i in range(self.n):
            for j in range(self.n):
                if self[i, j] == TileType.ELBOW_SE:
                    elbows.add((i, j))
        return elbows

    def all_blank_spots(self):
        blanks = set()
        for i in range(self.n):
            for j in range(self.n):
                if self[i, j] == TileType.BLANK:
                    blanks.add((i, j))
        return blanks

    def droop_moves(self):
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

    def do_droop_move(self, move):
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

    def pop_op(self):
        # --- STEP 0 --- #
        D = self.copy()
        # check if D has a blank tile (i.e., the coxeter length of D.w is zero)
        if self.perm.inv == 0:
            return D

        # find the first row r with a blank tile
        r = min([i for i in range(D.n) for j in range(D.n) if D[i, j] == TileType.BLANK])

        # initialize the mark X at the blank tile (x,y)
        x = r
        y = max([j for i in range(D.n) for j in range(D.n) if D[i, j] == TileType.BLANK and i == r])

        while True:
            # --- STEP 1 --- #
            # move the mark to the rightmost blank tile in the block
            j = y
            while j < D.n and D[x, j] == TileType.BLANK:
                j += 1
            y = j - 1

            # --- STEP 2 --- #
            # find first j-elbow (x_,y+1) with x_>x, if not found then set x_=0 (in which case p==y+1)
            x_ = 0
            for i in range(x + 1, D.n):
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
        for i in range(D.n - 1, x, -1):
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
        return D, (a + 1, r + 1)

    def resize(self, new_n):
        D = self.copy()
        D.grid.resize((new_n, new_n), refcheck=False)
        for i in range(D.n, new_n):
            for j in range(D.n, new_n):
                D.grid[i, j] = TileType.TBD
        D.n = new_n
        D.rebuild()
        return D

    def inverse_pop_op(self, a, r):
        D = self.copy()
        # check if D has a blank tile (i.e., the coxeter length of D.w is zero)
        if D.n <= a or D.n < r:
            new_n = max(D.n, a + 1, r)
            D = D.resize(new_n)
        # find first elbow in column a
        x_ = D.n - 1
        while x_ >= 0 and D[x_, a] != TileType.ELBOW_SE:
            x_ -= 1
        if x_ == -1:
            raise ValueError("No elbow found in specified column for inverse pop operation")
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
            elif D[z, y] == TileType.CROSS:
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
                elif D[z, y] == TileType.CROSS:
                    D.grid[z, y] = TileType.TBD
                    D.grid[z, y + 1] = TileType.CROSS
                elif D[z, y] == TileType.VERT and D[z, y + 1] == TileType.ELBOW_SE:
                    D.grid[z, y] = TileType.ELBOW_SE
                    D.grid[z, y + 1] = TileType.CROSS

            D.rebuild()

        D.rebuild()
        return D

        # r = min([i for i in range(D.n) for j in range(D.n) if D[i, j] == TileType.BLANK])

        # # initialize the mark X at the blank tile (x,y)
        # x = r
        # y = max([j for i in range(D.n) for j in range(D.n) if D[i, j] == TileType.BLANK and i == r])

        # while True:
        #     # --- STEP 1 --- #
        #     # move the mark to the rightmost blank tile in the block
        #     j = y
        #     while j < D.n and D[x, j] == TileType.BLANK:
        #         j += 1
        #     y = j - 1

        #     # --- STEP 2 --- #
        #     # find first j-elbow (x_,y+1) with x_>x, if not found then set x_=0 (in which case p==y+1)
        #     x_ = 0
        #     for i in range(x + 1, D.n):
        #         if D[i, y + 1] == TileType.ELBOW_NW:
        #             x_ = i
        #             break
        #     if x_ == 0:  # p==y+1
        #         break

        #     for z in range(x + 1, x_):
        #         if D[z, y] == TileType.BLANK:
        #             D.grid[z, y] = TileType.VERT
        #             D.grid[z, y + 1] = TileType.BLANK
        #         elif D[z, y] == TileType.CROSS:
        #             continue
        #         elif D[z, y] == TileType.HORIZ:
        #             D.grid[z, y] = TileType.CROSS
        #             D.grid[z, y + 1] = TileType.HORIZ
        #         elif D[z, y] == TileType.VERT:
        #             continue
        #         elif D[z, y] == TileType.ELBOW_SE:
        #             D.grid[z, y] = TileType.VERT
        #             D.grid[z, y + 1] = TileType.ELBOW_SE
        #         elif D[z, y] == TileType.ELBOW_NW:
        #             D.grid[z, y] = TileType.CROSS
        #             D.grid[z, y + 1] = TileType.ELBOW_NW
        #     D.grid[x, y] = TileType.ELBOW_SE  # NW-corner
        #     D.grid[x_, y + 1] = TileType.BLANK  # SE-corner
        #     if D[x_, y] == TileType.ELBOW_SE:  # SW-corner
        #         D.grid[x_, y] = TileType.VERT
        #     else:  # D.grid[x_,y] == TileType.HORIZ
        #         D.grid[x_, y] = TileType.ELBOW_NW
        #     if D[x, y + 1] == TileType.ELBOW_SE:  # NE-corner
        #         D.grid[x, y + 1] = TileType.HORIZ
        #     else:  # D.grid[x,y+1] == TileType.VERT
        #         D.grid[x, y + 1] = TileType.ELBOW_NW

        #     # move the mark X to the SE-corner of U
        #     x = x_
        #     y = y + 1

        # # --- STEP 3 --- #
        # a = y  # where (x,y) is the final position of the mark X

        # x_ = 0
        # for i in range(D.n - 1, x, -1):
        #     if D[i, y] == TileType.ELBOW_SE:
        #         x_ = i
        #         break

        # # copied from above
        # for z in range(x + 1, x_):
        #     if D[z, y] == TileType.BLANK:
        #         D.grid[z, y] = TileType.VERT
        #         D.grid[z, y + 1] = TileType.BLANK
        #     elif D[z, y] == TileType.CROSS:
        #         continue
        #     elif D[z, y] == TileType.HORIZ:
        #         D.grid[z, y] = TileType.CROSS
        #         D.grid[z, y + 1] = TileType.HORIZ
        #     elif D[z, y] == TileType.VERT:
        #         continue
        #     elif D[z, y] == TileType.ELBOW_SE:
        #         D.grid[z, y] = TileType.VERT
        #         D.grid[z, y + 1] = TileType.ELBOW_SE
        #     elif D[z, y] == TileType.ELBOW_NW:
        #         D.grid[z, y] = TileType.CROSS
        #         D.grid[z, y + 1] = TileType.ELBOW_NW

        # D.grid[x, y] = TileType.ELBOW_SE  # NW-corner
        # D.grid[x_, y + 1] = TileType.ELBOW_SE  # SE-corner
        # D.grid[x_, y] = TileType.VERT  # SW-corner
        # if D[x, y + 1] == TileType.ELBOW_SE:  # NE-corner
        #     D.grid[x, y + 1] = TileType.HORIZ
        # else:  # D.grid[x,y+1] == TileType.VERT
        #     D.grid[x, y + 1] = TileType.ELBOW_NW

        # D.rebuild()
        # return D, (a + 1, r + 1)

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
            work_bpd, (reflection, row) = work_bpd.pop_op()
            rows[row - 1] = rows[row - 1] + [reflection]

        return RCGraph([tuple(r) for r in rows])

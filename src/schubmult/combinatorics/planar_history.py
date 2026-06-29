from functools import cached_property

import numpy as np


class Tile:
    NORTH = 0
    EAST = 1
    SOUTH = 2
    WEST = 3

    def __init__(self, edges: set[tuple[int, int]] | None = None):
        self.edges = edges if edges is not None else set()

    def __hash__(self):
        return hash(frozenset(self.edges))

    def __eq__(self, other):
        return self.edges == other.edges

class PlanarHistory:

    CROSS = type("Tile", (Tile,), {"__str__": lambda self: "┼", "__repr__": lambda self: "┼"})(edges={(Tile.SOUTH, Tile.NORTH), (Tile.WEST, Tile.EAST)})
    BUMP = type("Tile", (Tile,), {"__str__": lambda self: "*", "__repr__": lambda self: "*"})(edges={(Tile.SOUTH, Tile.EAST), (Tile.WEST, Tile.NORTH)})
    EMPTY = type("Tile", (Tile,), {"__str__": lambda self: " ", "__repr__": lambda self: " "})(edges=set())

    def __init__(self, grid: np.ndarray):
        self._grid = grid.copy()
        self._rows = grid.shape[0]
        self._cols = grid.shape[1]

    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    def __getitem__(self, key):
        return self._grid[key]

    @cached_property
    def grid(self):
        return self._grid.copy()

    @property
    def _ne_grid(self):
        return self._grid

    def _tile_ne_contribution(self, tile: Tile) -> int:
        """Contribution of a tile to the NE pipe count recurrence."""
        if tile == self.EMPTY:
            return 0
        if tile == self.CROSS:
            return 2
        return 1

    def _perm_word_letter(self, pipes_northeast: int) -> int:
        """Convert a NE pipe count to a simple-reflection letter."""
        return int(pipes_northeast - 1)

    @property
    def perm_word(self):
        cross = self.CROSS
        diff = np.zeros((self.rows, self.cols), dtype=int)
        for i in range(self.rows):
            for j in range(self.cols):
                diff[i, j] = self._tile_ne_contribution(self._ne_grid[i, j])
        # Create r array with shape (self.rows+1, self.cols+1)
        r = np.zeros((self.rows + 1, self.cols + 1), dtype=int)

        for i in range(1, self.rows + 1):
            for j in range(1, self.cols + 1):
                r[i, j] = r[i - 1, j - 1] + diff[i - 1, j - 1]

        # Pre-compute all cross positions and their pipes_northeast values
        cross_positions = np.argwhere(self._ne_grid == cross)
        word = []
        if len(cross_positions) > 0:
            # Sort by column first, then by row descending (for correct swap order)
            sort_indices = np.lexsort((-cross_positions[:, 0], cross_positions[:, 1]))
            cross_positions = cross_positions[sort_indices]
            # Get pipes_northeast for each cross position
            pipes_northeast_values = r[cross_positions[:, 0] + 1, cross_positions[:, 1] + 1]
            # Apply swaps sequentially
            for pipes_northeast in pipes_northeast_values:
                word = [*word, self._perm_word_letter(int(pipes_northeast))]

        return tuple(word)

import numpy as np


class Tile:
    NORTH = 0
    EAST = 1
    SOUTH = 2
    WEST = 3

    def __init__(self, edges: set[tuple[int, int]]):
        self.edges = edges

    def __eq__(self, other):
        if not isinstance(other, Tile):
            return NotImplemented
        return self.edges == other.edges

class PlanarHistory:
    def __init__(self, grid: np.ndarray):
        self._grid = grid.copy()
        self.rows = grid.shape[0]
        self.cols = grid.shape[1]

    @property
    def perm_word(self):
        cross = Tile(edges={(Tile.SOUTH, Tile.NORTH), (Tile.WEST, Tile.EAST)})
        diff = np.ones((self.rows, self.cols), dtype=int)
        diff[self._grid == Tile(edges=set())] = 0
        diff[self._grid == cross] = 2
        # Create r array with shape (self.rows+1, self.cols+1)
        r = np.zeros((self.rows + 1, self.cols + 1), dtype=int)

        for i in range(1, self.rows + 1):
            for j in range(1, self.cols + 1):
                r[i, j] = r[i - 1, j - 1] + diff[i - 1, j - 1]

        # Pre-compute all cross positions and their pipes_northeast values
        cross_positions = np.argwhere(self._grid == cross)
        word = []
        if len(cross_positions) > 0:
            # Sort by column first, then by row descending (for correct swap order)
            sort_indices = np.lexsort((-cross_positions[:, 0], cross_positions[:, 1]))
            cross_positions = cross_positions[sort_indices]
            # Get pipes_northeast for each cross position
            pipes_northeast_values = r[cross_positions[:, 0] + 1, cross_positions[:, 1] + 1]
            # Apply swaps sequentially
            for pipes_northeast in pipes_northeast_values:
                # THIS LINE IS WRONG
                # if small_perm[pipes_northeast - 2] < small_perm[pipes_northeast - 1]:
                word = [*word, int(pipes_northeast - 1)]

        return tuple(word)

from functools import cached_property

import numpy as np

from schubmult.combinatorics.permutation import Permutation
from schubmult.utils._grid_print import GridPrint

from .planar_history import PlanarHistory


class PipeDream(PlanarHistory, GridPrint):
    # def __str__(self) -> str:
    #     symbols = {Tile(edges=set()): " ", PipeDream.CROSS: "┼", PipeDream.BUMP: "*"}
    #     return symbols.get(self, "?")

    # def __repr__(self) -> str:
    #     symbols = {Tile(edges=set()): " ", PipeDream.CROSS: "┼", PipeDream.BUMP: "*"}
    #     return symbols.get(self, "?")

    def __init__(self, grid: np.ndarray):
        super().__init__(grid)

    @cached_property
    def _native_cross_word(self) -> tuple[int, ...]:
        """Pipe-dream word in native left->top orientation.

        Read crosses row-by-row from top to bottom, and right-to-left inside each row.
        A cross at zero-based position ``(r, c)`` contributes the letter ``r + c + 1``.
        """
        cross_positions = np.argwhere(self._grid == self.CROSS)
        if len(cross_positions) == 0:
            return ()
        sort_indices = np.lexsort((-cross_positions[:, 1], cross_positions[:, 0]))
        cross_positions = cross_positions[sort_indices]
        return tuple(int(r + c + 1) for r, c in cross_positions)

    @property
    def perm_word(self):
        return self._native_cross_word

    @property
    def is_reduced(self):
        return self.perm.inv == len(self.perm_word)

    @property
    def perm(self):
        # Use native pipe-dream orientation (left -> top) for permutation extraction.
        return Permutation.hecke_ref_product(*self._native_cross_word)

    def to_rc_graph(self):
        from schubmult import RCGraph

        rows = []
        for i in range(self.rows):
            row = []
            for j in range(self.cols - i):
                if self[i, j] == self.CROSS:
                    row.append(i + j + 1)
            rows.append(tuple(reversed(row)))
        return RCGraph(rows)

    def to_wc_graph(self):
        from schubmult import WCGraph

        rows = []
        for i in range(self.rows):
            row = []
            for j in range(self.cols - i):
                if self[i, j] == self.CROSS:
                    row.append(i + j + 1)
            rows.append(tuple(reversed(row)))
        return WCGraph(rows)

    @classmethod
    def from_rc_graph(cls, rc_graph):
        grid = np.full((len(rc_graph.perm), len(rc_graph.perm)), cls.EMPTY, dtype=object)
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1] - i):
                if rc_graph.has_element(i + 1, j + 1):
                    grid[i, j] = cls.CROSS
                else:
                    grid[i, j] = cls.BUMP
        return cls(grid)

    @classmethod
    def from_wc_graph(cls, wc_graph):
        grid = np.full((len(wc_graph.perm), len(wc_graph.perm)), cls.EMPTY, dtype=object)
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1] - i):
                if wc_graph.has_element(i + 1, j + 1):
                    grid[i, j] = cls.CROSS
                else:
                    grid[i, j] = cls.BUMP
        return cls(grid)

    @property
    def _ne_grid(self):
        # Pipe dreams are stored with pipes entering from the left and exiting at the top.
        # Rotate to convert that boundary orientation to NE-planar-history convention
        # (enter bottom, exit right) before applying PlanarHistory word extraction.
        return np.rot90(self._grid, 2)

    def co_pipe_dream(self):
        new_grid = self.grid.copy()
        new_grid[:] = self.EMPTY
        for i in range(1, self.grid.shape[0] + 1):
            for j in range(1, self.grid.shape[1] + 1 - i):
                if self[i - 1, j - 1] == self.CROSS:
                    new_grid[self.rows - 1 - (i + j - 1), j - 1] = self.BUMP
                if self[i - 1, j - 1] == self.BUMP:
                    new_grid[self.rows - 1 - (i + j - 1), j - 1] = self.CROSS
                # new_grid[i + j - 1, j - 1] = self.CROSS if self[i - 1, j - 1] == self.BUMP else (self.BUMP if self[i - 1, j - 1] == self.CROSS else self.EMPTY)
        return PipeDream(new_grid)

    def co_stinkbat_pipe_dream(self):
        new_grid = self.grid.copy()
        new_grid[:] = self.EMPTY
        for i in range(1, self.grid.shape[0] + 1):
            for j in range(1, self.grid.shape[1] + 1 - i):
                if self[i - 1, j - 1] == self.CROSS:
                    new_grid[self.rows - 1 - (i + j - 1), j - 1] = self.CROSS
                if self[i - 1, j - 1] == self.BUMP:
                    new_grid[self.rows - 1 - (i + j - 1), j - 1] = self.BUMP
                # new_grid[i + j - 1, j - 1] = self.CROSS if self[i - 1, j - 1] == self.BUMP else (self.BUMP if self[i - 1, j - 1] == self.CROSS else self.EMPTY)
        return PipeDream(new_grid)

    def inverse_pipe_dream(self):
        new_grid = self.grid.copy()
        new_grid[:] = self.EMPTY
        for i in range(1, self.grid.shape[0] + 1):
            for j in range(1, self.grid.shape[1] + 1 - i):
                if self[self.rows - 1 - (i + j - 1), j - 1] == self.BUMP:
                    new_grid[i - 1, j - 1] = self.CROSS
                else:
                    new_grid[i - 1, j - 1] = self.BUMP
                # new_grid[i + j - 1, j - 1] = self.CROSS if self[i - 1, j - 1] == self.BUMP else (self.BUMP if self[i - 1, j - 1] == self.CROSS else self.EMPTY)
        return PipeDream(new_grid)

    def reflect_vertically(self):
        new_grid = np.flipud(self.grid.copy())
        return PipeDream(new_grid)

    def __hash__(self):
        return hash(tuple(map(tuple, self.grid)))

    def __eq__(self, other):
        if not isinstance(other, PipeDream):
            return NotImplemented
        return np.array_equal(self._grid, other._grid)

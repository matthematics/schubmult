from functools import cached_property
from typing import Any, Optional

import numpy as np

from schubmult.utils._grid_print import GridPrint

from .crystal_graph import CrystalGraph
from .nilplactic import NilPlactic
from .perm_lib import Permutation
from .plactic import Plactic
from .rc_graph import RCGraph

# specific skew tableaux
# subword
# dual
# recti/fies to a specific subword
#subword subtableaux

# Dominant formula: number of standard tableaux of shape w/mu that recitfy to a subword tableau of shape u for a fixed tableau of shape w
# skew tableau behave dominant correctly
# w/mu subword tableau that are a specific standard tableaux for v
class RootTableau(CrystalGraph, GridPrint):
    """
    Root tableau with dual knuth equivalence
    """

    def __hash__(self) -> int:
        return hash(self._hasher)


    @classmethod
    def root_insert_rsk(cls, reduced_word, compatible_seq):
        _perm = Permutation.ref_product(*reduced_word)
        word, word2 = (), ()
        spunkle = len(_perm)
        w0 = Permutation.w0(spunkle)
        rev_word = [w0[a-1] for a in reduced_word]

        for idx, letter in enumerate(rev_word):
            letter2 = idx + 1
            word, word2 = NilPlactic._ed_insert_rsk(word, word2, int(letter), int(letter2) if letter2 is not None else None)
        num_rows = len(word2)
        num_cols = max(len(r) for r in word2)
        grid = np.empty((num_rows, num_cols), dtype=object)
        for r in range(num_rows):
            for c in range(num_cols):
                if c < len(word2[r]):
                    # store pair (root_tuple, letter) as before
                    grid[r, c] = (_perm.right_root_at(word2[r][c] - 1, word=reduced_word), compatible_seq[word2[r][c] - 1])
                else:
                    grid[r, c] = None

        return cls(grid)

    # skew tableaux are subword
    @classmethod
    def from_rc_graph(cls, rc: RCGraph):
        reduced_word = rc.perm_word
        compatible_seq = []
        for i in range(len(rc)):
            compatible_seq.extend([i+1] * len(rc[i]))
        return cls.root_insert_rsk(reduced_word, compatible_seq)

    def delete_box(self, i, j):
        pass

    def __getitem__(self, key: Any) -> Any:
        return self._root_grid[key]

    @cached_property
    def rows(self):
        return self._root_grid.shape[0]

    @cached_property
    def cols(self):
        return self._root_grid.shape[1]

    # @property
    # def compatible_sequence(self):
    #     return self._plactic.reverse_rsk(self._index_tableau)

    @property
    def is_valid(self):
        return self.rc_graph.is_valid

    @property
    def reduced_word(self):
        def root_shift(d):
            return np.vectorize(lambda x: (Permutation.ref_product(d).act_root(*x[0]), x[1]) if x is not None else None, otypes=[object])
        grid = self._root_grid.copy()
        word = []
        def _recurse_grid():
            nonlocal grid, word
            max_r, max_c = grid.shape
            for r in range(max_r - 1, -1, -1):
                if grid[r, 0] is not None:
                    for c in range(max_c - 1, -1, -1):
                        if grid[r, c] is not None:
                            d = grid[r, c][0][0]
                            word.append(d)
                            grid[r, c] = None
                            grid = root_shift(d)(grid)
                            return True
            return False
        while _recurse_grid():
            pass
        return tuple(reversed(word))
    # @property
    # def reduced_word(self):
    #     return self._red_plactic.reverse_rsk(self._index_tableau)

    def __init__(self, grid=None):
        self._root_grid = grid.copy()
        self._hasher = tuple(tuple(tuple(b) for b in a if b is not None) for a in self._root_grid if a is not None)


    @property
    def weight_tableau(self):
        _word = tuple([a for a in row if a is not None] for row in self._root_grid)
        return Plactic(_word)

    def epsilon(self, index):
        return self._weight_tableau.epsilon(index)

    def raising_operator(self, index):
        new_plactic = self._weight_tableau.raising_operator(index)
        if new_plactic is None:
            return None
        return RootTableau(self._base_word, new_plactic)

    def lowering_operator(self, index):
        new_plactic = self._weight_tableau.lowering_operator(index)
        if new_plactic is None:
            return None
        new_tab = RootTableau(self._base_word, new_plactic)
        if new_tab.rc_graph is None or not new_tab.rc_graph.is_valid:
            return None
        return new_tab

    @property
    def rc_graph(self):
        rows = []
        for i in range(len(self._base_word)):
            if i == 0 or (i > 0 and self._base_word[i] > self._base_word[i - 1]):
                rows.append([])
            rows[-1].append(self._base_word[i])
        raise_seq = self._weight_tableau.to_highest_weight(length=self.crystal_length())[1]
        try:
            rc = RCGraph([tuple(row) for row in rows]).normalize()
            return rc.reverse_raise_seq(raise_seq)
        except ValueError:
            return None

    def crystal_length(self):
        return len(self._perm.trimcode)

    @property
    def crystal_weight(self):
        return self._plactic.crystal_weight()

    # def reduced_word(self):
    #     recording_tableau = self.index_tableau
    #     permutation_of_roots = self.reverse_rsk(recording_tableau)
    #     roots = list(reversed([self._perm.right_root_at(i - 1) for i in permutation_of_roots]))
    #     word = []
    #     for i in range(len(roots)):
    #         word.append(roots[0][0])
    #         sref = Permutation.ref_product(roots[0][0])
    #         roots = [sref.act_root(*r) for r in roots[1:]]
    #     return tuple(reversed(word))




    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RootTableau):
            return False
        return self._hasher == other._hasher

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
# rectifies to a specific subword
#subword subtableaux

# Dominant formula: number of standard tableaux of shape w/mu that recitfy to a subword tableau of shape u for a fixed tableau of shape w
# skew tableau behave dominant correctly
# w/mu subword tableau that are a specific standard tableaux for v
class RootTableau(CrystalGraph, GridPrint):
    """
    Root tableau with dual knuth equivalence
    """

    def __hash__(self) -> int:
        return hash(tuple(self._base_word)) ^ hash(self._weight_tableau)

    # edelman greene highest weight
    @classmethod
    def root_rsk_insert(cls, reduced_word, compatible_seq):
        grid = np.empty((len(reduced_word), len(reduced_word)), dtype=object)
        def _recurse(index):
            nonlocal grid
            if index == len(reduced_word):
                return cls(
                    grid,
                )
            if index == 0:
                d = reduced_word[0]
                grid[0, 0] = ((d, d + 1), compatible_seq[0])
                return _recurse(index + 1)
            d = reduced_word[index]
            root_shift = np.vectorize(lambda x: (Permutation.ref_product(d).act_root(*x[0]), x[1]) if x is not None else None, otypes=[object])
            grid = root_shift(grid)
            to_insert = compatible_seq[index]
            root = (d, d +  1)
            def _root_compare(root1, root2):
                if root1[0] == root2[0]:
                    return 1
                if root1[1] == root2[1]:
                    return 1
                if root1[0] == root2[1] or root1[1] == root2[0]:
                    return -1
                return 0

            for i in range(len(reduced_word)):
                last_box = max((j for j in range(len(reduced_word)) if grid[i, j] is not None), default=-1)
                if last_box == -1:
                    grid[i, 0] = (root, to_insert)
                    return _recurse(index + 1)
                good = True
                for k in range(last_box + 1):
                    (a, b), entry = grid[i, k]
                    if to_insert < entry:
                        # bump here
                        grid[i, k] = (root, to_insert)
                        (root, to_insert) = (a, b), entry
                        break
                    if _root_compare(root, (a, b)) == -1:
                        good = False
                        break
                if not good:
                    continue
                grid[i, last_box + 1] = (root, to_insert)
                return _recurse(index + 1)
            grid[i, 0] = (root, to_insert)
            return _recurse(index + 1)

        return _recurse(0)
            # how do we insert the simple reflection


    # skew tableaux are subword
    @classmethod
    def from_rc_graph(cls, rc: RCGraph):
        rc_hw, raise_seq = rc.to_highest_weight()
        weight_tableau = rc_hw.weight_tableau
        base_word = rc.perm_word
        return cls(base_word, weight_tableau.reverse_raise_seq(raise_seq))

    def delete_box(self, i, j):
        if j >= self.weight_tableau.shape[i]:
            return None
        index = sum(self.weight_tableau.shape[:i]) + j
        
        new_word = [list(row) for row in self.weight_tableau._word]
        new_word[i][j] = None
        new_w_tab = Plactic(new_word).up_jdt_slide(i, j)
        base_word = self._base_word[:index] + self._base_word[index + 1 :]
        return RootTableau(base_word, new_w_tab)

    def __getitem__(self, key: Any) -> Any:
        return self._root_grid[key]
        # if isinstance(key, tuple):
        #     i, j = key
        #     if i >= self.weight_tableau.rows:
        #         return self.weight_tableau[i - self.weight_tableau.rows, j]
        #     if j >= self.weight_tableau.shape[i]:
        #         return None
        #     index = sum(self.weight_tableau.shape[:i]) + j
        #     try:
        #         return self.perm.right_root_at(index, word=self.base_word)
        #     except IndexError:
        #         return None
        # is_slice = isinstance(key, slice)
        # if is_slice:
        #     start = key.start or 0
        #     stop = key.stop or len(self.base_word)
        #     step = key.step or 1
        #     roots = []
        #     for idx in range(start, stop, step):
        #         roots.append(self.perm.right_root_at(idx, word=self.base_word))
        #     return tuple(roots)
        # raise TypeError(f"Invalid key type {type(key)}")

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
        # check if the tableau is a valid root tableau
        # i.e. for each box (r,c) with entry (a,b), we have a <= b
        return self.rc_graph.is_valid

    @property
    def perm(self):
        return self._perm

    @property
    def weight_tableau(self):
        return self._weight_tableau

    @property
    def base_word(self):
        return self._base_word

    # @property
    # def reduced_word(self):
    #     return self._red_plactic.reverse_rsk(self._index_tableau)

    def __init__(self, grid=None):
        self._root_grid = grid


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
        return self._word == other._word
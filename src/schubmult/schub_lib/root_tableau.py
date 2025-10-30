
from functools import cached_property
from typing import Any, Optional

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

    @classmethod
    def from_rc_graph(cls, rc: RCGraph):
        rc_hw, raise_seq = rc.to_highest_weight()
        weight_tableau = rc_hw.weight_tableau
        base_word = rc_hw.perm_word
        return cls(base_word, weight_tableau.reverse_raise_seq(raise_seq))

    def delete_box(self, index):
        i = index - 1
        wd = [*self._red_plactic]
        root = self.perm.right_root_at(i, word=wd)
        wd.pop(i)
        pm = Permutation.ref_product(*wd)
        if pm.inv != len(wd):
            return None
        return RootTableau(reduced_word=wd), root

    def __getitem__(self, key: Any) -> Any:
        return self._weight_tableau[key]

    @cached_property
    def rows(self):
        return self.weight_tableau.rows + 1

    @cached_property
    def cols(self):
        return max(self.weight_tableau.cols,len(self.base_word))

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

    def __init__(self, base_word, weight_tableau):
        self._weight_tableau = weight_tableau
        self._base_word = base_word
        self._perm = Permutation.ref_product(*self._base_word)


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
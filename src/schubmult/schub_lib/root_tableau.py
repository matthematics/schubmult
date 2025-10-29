
from functools import cached_property
from typing import Any, Optional

from .crystal_graph import CrystalGraphDual
from .nilplactic import NilPlactic
from .perm_lib import Permutation
from .plactic import Plactic


# dual
class RootTableau(Plactic):
    """
    Root tableau with dual knuth equivalence
    """

    def __getitem__(self, key: Any) -> Any:
        value = super().__getitem__(key)
        if isinstance(value, int):
            return self.perm.right_root_at(value)
        return tuple([self.perm.right_root_at(i) for i in value])


    def __init__(self, reduced_word=()):
        perm = Permutation.ref_product(*reduced_word)
        root_seq = tuple(perm.right_root_at(i-1, word=list(reversed(reduced_word))) for i in range(1, len(perm)))
        permutation_of_roots = [root_seq.index(perm.right_root_at(i) for i in range(1, len(perm))) + 1]
        plactic = Plactic().rs_insert(*permutation_of_roots)
        super().__init__(plactic._word)
        self._perm = perm

    @property
    def perm(self) -> Permutation:
        return self._perm

    # @property
    # def perm_tab_word(self):
    #     return tuple(reversed(self._nilplactic.row_word))

    # @property
    # def perm(self):
    #     return ~(self._nilplactic.perm)

    # def _root(self, i):
    #     return self.perm.right_root_at(i-1, word=self.perm_tab_word)

    # # 1-indexed
    # #1515 Zurika
    # def root_label(self, i):
    #     return self.row_word[len(self.row_word) - i]

    # @classmethod
    # def from_rc(cls, reduced_word, compat_seq):
    #     _perm = ~Permutation.ref_product(*reduced_word)
    #     _root_seq = tuple([_perm.right_root_at(i-1, word=list(reversed(reduced_word))) for i in range(1, len(_perm))])
    #     _compat_seq = tuple(compat_seq)
    #     _hw_seq = []]
    #     index = 0
    #     last_min_desc = 100
    #     cur_value = 0
    #     working_perm = _perm
    #     red_word = list(reduced_word)
    #     while index < len(_root_seq):
    #         d = min(~working_perm).descents
            
        

    # def crystal_length(self) -> int:
    #     return len(self.perm.trimcode)
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RootTableau):
            return False
        return self._plactic == other._plactic and self._nilplactic == other._nilplactic
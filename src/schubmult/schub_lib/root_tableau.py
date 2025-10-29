
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
            return self.perm.right_root_at(value - 1)
        return tuple([self.perm.right_root_at(i - 1) for i in value])


    def __init__(self, reduced_word=()):
        perm = Permutation.ref_product(*reduced_word)
        assert perm.inv == len(reduced_word), "Reduced word does not match permutation length"
        root_seq = tuple(perm.right_root_at(i, word=reduced_word) for i in range(perm.inv))
        permutation_of_roots = [root_seq.index(perm.right_root_at(i, word=reduced_word)) + 1 for i in range(perm.inv)]
        plactic = Plactic().rs_insert(*permutation_of_roots)
        super().__init__(plactic._word)
        self._perm = perm

    @property
    def perm(self) -> Permutation:
        return self._perm

    
    
    def reduced_word(self, recording_tableau):
        permutation_of_roots = self.reverse_rsk(recording_tableau)
        roots = list(reversed([self._perm.right_root_at(i - 1) for i in permutation_of_roots]))
        word = []
        for i in range(len(roots)):
            word.append(roots[0][0])
            sref = Permutation.ref_product(roots[0][0])
            roots = [sref.act_root(*r) for r in roots[1:]]
        return tuple(reversed(word))



    def is_compatible_sequence(self, seq):
        if any(seq[i-1] > seq[i] for i in range(1, len(seq))):
            return False
        insert_seq = [seq[self.row_word[i] - 1] for i in range(len(self.row_word))]
        plactic = Plactic().rs_insert(*insert_seq)
        rw = self.reduced_word(Plactic.superstandard(self.shape))
        test_word = plactic.reverse_rsk(Plactic.superstandard(plactic.shape))
        for letter1, letter2 in zip(test_word, rw):
            if letter1 > letter2:
                return False
        return True

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RootTableau):
            return False
        return self._plactic == other._plactic and self._nilplactic == other._nilplactic
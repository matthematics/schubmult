from functools import cache

from schubmult.schub_lib.permutation import Permutation, uncode
from schubmult.symbolic import Symbol, sstr

from ..schubert_ring import Sx
from .free_algebra_basis import FreeAlgebraBasis


class SchubertSchurBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return len(x) == 2 and isinstance(x[0], list | tuple) and isinstance(x[1], Permutation | list | tuple)

    @classmethod
    def as_key(cls, x):
        return (tuple(x[0]), Permutation(x[1]))

    zero_monom = ((), Permutation([]))

    @classmethod
    @cache
    def coproduct(cls, key):
        from ...utils._mul_utils import _tensor_product_of_dicts_first
        from .schubert_basis import SchubertBasis

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(SchubertBasis.transition(cls)(x[0]), SchubertBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda x: SchubertBasis.coproduct(x), cls.transition_schubert(*key)),
        )

    @classmethod
    def transition_schubert(cls, lambd, perm):
        from ..schubert_ring import SingleSchubertRing
        from ..variables import MaskedGeneratingSet

        if lambd[-1] == 0:
            return {(perm, len(lambd)): 1}
        numvars = len(lambd)
        extra = lambd[-1] + len(lambd) - 1
        dom = uncode([numvars] * extra)
        grass_perm = uncode(lambd) * dom
        w0 = uncode(list(range(numvars - 1, 0, -1)))
        lower_perm = perm * w0
        dom_perm = uncode(([numvars] * extra) + list(range(numvars - 1, 0, -1)))
        shifted_ring = SingleSchubertRing(MaskedGeneratingSet(Sx([]).ring.genset, list(range(1, extra + 1))))
        start_schub = Sx(grass_perm)
        start_schub *= shifted_ring(lower_perm).in_SEM_basis()
        return {(k * ~dom_perm, numvars): v for k, v in start_schub.items()}

    @classmethod
    def transition_word(cls, lambd, perm):
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis), cls.transition_schubert(lambd, perm))

    @classmethod
    def transition(cls, other_basis):
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: x
        return None

    @classmethod
    def printing_term(cls, k):
        return Symbol(f"SS{sstr(k)}")

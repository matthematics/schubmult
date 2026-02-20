from itertools import zip_longest

from schubmult.schub_lib.permutation import Permutation, uncode
from schubmult.symbolic import S, Symbol, sstr
from schubmult.utils.perm_utils import add_perm_dict

from ..schubert.schubert_ring import Sx
from .free_algebra_basis import FreeAlgebraBasis


class _SeparatedDescentsBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return len(x) == 3 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], Permutation | list | tuple) and isinstance(x[2], int) and (len(x[1]) <= x[2])

    @classmethod
    def as_key(cls, x):
        return (Permutation(x[0]), Permutation(x[1]), x[2])

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        from .schubert_basis import SchubertBasis

        left = cls.transition_schubert(*key1)
        right = cls.transition_schubert(*key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(SchubertBasis.transition(cls), SchubertBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret

    @classmethod
    def transition_schubert(cls, perm0, perm1, numvars):
        from schubmult.symbolic.poly.variables import MaskedGeneratingSet

        from ..schubert.schubert_ring import SingleSchubertRing

        mu1_code = perm0.mul_dominant().trimcode
        mu2_code = list(range(cls.k - 1, 0, -1))
        mu_code = [a + b for a, b in zip_longest(mu1_code, mu2_code, fillvalue=0)]

        mu = uncode(mu_code)
        mu1 = uncode(mu1_code)
        mu2 = uncode(mu2_code)

        mul1 = perm0 * (~mu1)
        mul2 = perm1 * (~mu2)
        shifted_ring = SingleSchubertRing(MaskedGeneratingSet(Sx([]).ring.genset, list(range(1, len((~mu).trimcode) - cls.k + 2))))
        start_schub = Sx(mul1)
        start_schub *= shifted_ring(mul2).in_SEM_basis()
        return {(k1 * mu, numvars): v for k1, v in start_schub.items()}

    @classmethod
    def transition_word(cls, perm0, perm1, n):
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis), cls.transition_schubert(perm0, perm1, n))

    @classmethod
    def transition(cls, other_basis):
        from .schubert_basis import SchubertBasis
        from .schubert_schur_basis import SchubertSchurBasis
        from .word_basis import WordBasis

        if other_basis == cls:
            return lambda x: {x: S.One}
        if isinstance(other_basis, type) and other_basis.__name__ == cls.__name__ and other_basis.k != cls.k:
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_separated_descents(cls.k, *y), cls.transition_schubert(*x))
        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_schubert_schur(*y), cls.transition_schubert(*x))
        return None

    @classmethod
    def printing_term(cls, k):
        return Symbol(f"SepDesc{cls.k}{sstr(k)}")


def SeparatedDescentsBasis(k):
    return type("_SeparatedDescentsBasis", (_SeparatedDescentsBasis,), {"k": k, "zero_monom": (Permutation(()), Permutation([]), 0)})

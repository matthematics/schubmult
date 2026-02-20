import schubmult.rings.free_algebra as fa
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class ZBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        return tuple(x)

    @staticmethod
    def from_perm(perm, n):
        cd = perm.code
        if len(cd) < n:
            return None
        if 0 in cd[:n]:
            return None
        return cd[:n]

    @staticmethod
    def pare_schubert(perm):
        cd = [*perm.code]
        while len(cd) > 0 and cd[-1] == 0:
            cd = cd[:-1]
        if 0 in cd:
            return None
        return tuple(cd)

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        from .schubert_basis import SchubertBasis

        key11 = uncode([a - 1 for a in key1])
        key22 = uncode([a - 1 for a in key2])

        dct = SchubertBasis.product((key11, len(key1)), (key22, len(key2)), coeff)

        JB = fa.FreeAlgebra(basis=ZBasis)
        dct_out = JB.zero

        for (perm, n), v in dct.items():
            cd = [*perm.code]
            tup = [1] * n
            for i, a in enumerate(cd[:n]):
                tup[i] = a + 1
            dct_out[tuple(tup)] = v
        return dct_out

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "Z")

    @classmethod
    def transition(cls, other_basis):
        from .word_basis import WordBasis

        def trans(x):
            from schubmult.rings import ASx

            dct = ASx(uncode(x)).change_basis(WordBasis)
            dct_out = {}
            for k, v in dct.items():
                if 0 in k:
                    continue
                k2 = k
                dct_out[k2] = dct_out.get(k2, S.Zero) + v
            return dct_out

        if other_basis == ZBasis:
            return lambda x: {x: S.One}
        if other_basis == WordBasis:
            return trans

        return lambda x: FreeAlgebraBasis.compose_transition(lambda k: WordBasis.transition(other_basis)(k), trans(x))

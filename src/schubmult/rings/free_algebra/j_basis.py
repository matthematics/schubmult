import schubmult.rings.free_algebra as fa
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class JBasis(FreeAlgebraBasis):
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

    @classmethod
    def coproduct(cls, key):
        return cls.bcoproduct(key)

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "J")

    @classmethod
    def transition(cls, other_basis):
        from .schubert_basis import SchubertBasis
        from .word_basis import WordBasis

        def trans(x):
            t = S.One
            ASx = fa.FreeAlgebra(basis=SchubertBasis)
            if 0 in x:
                raise ValueError("Transition from JBasis to other basis is only implemented for keys with no zeros.")
            dct = ASx(uncode(x)).change_basis(WordBasis)
            dct_out = {}
            for k, v in dct.items():
                k2 = tuple(a for a in k if a != 0)
                extra_coeff = t ** (len(k) - len(k2)) if len(k) - len(k2) > 0 else S.One
                dct_out[k2] = dct_out.get(k2, S.Zero) + v * extra_coeff
            return dct_out

        if other_basis == JBasis:
            return lambda x: {x: S.One}
        if other_basis == WordBasis:
            return trans

        return lambda x: FreeAlgebraBasis.compose_transition(lambda k: WordBasis.transition(other_basis)(k), trans(x))

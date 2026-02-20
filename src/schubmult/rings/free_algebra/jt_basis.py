from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S, Symbol, sstr, sympy_Mul

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class JTBasis(FreeAlgebraBasis):
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

    @staticmethod
    def normalize_dct(dct):
        dct_out = {}
        for k, v in dct.items():
            k2 = tuple([a for a in k if a != 0])
            pw = len(k) - len(k2)
            k2 = tuple([0] * pw + list(k2))
            if k2 not in dct_out:
                dct_out[k2] = S.Zero
            dct_out[k2] += v
        return dct_out

    zero_monom = ()
    t = Symbol("t")

    @classmethod
    def printing_term(cls, k):
        return sympy_Mul(JTBasis.t ** k[1], GenericPrintingTerm((k[0], "t"), "JT"))

    @classmethod
    def transition(cls, other_basis):
        from .word_basis import WordBasis

        @cache
        def trans(x2):
            from schubmult import ASx

            if len(x2) != 2:
                raise ValueError("JTBasis transition expects a tuple of length 2, got: " + sstr(x2))
            x = x2[0]
            n = x2[1]  # noqa: F841
            return JTBasis.normalize_dct(ASx(uncode(list(x)), len(x)).change_basis(WordBasis))

        if other_basis == JTBasis:
            return lambda x: {x: S.One}
        if other_basis == WordBasis:
            return trans
        return lambda x: FreeAlgebraBasis.compose_transition(lambda k: WordBasis.transition(other_basis)(k), trans(x))

from functools import cache

from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S, Symbol, sstr, sympy_Mul

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class JTBasis(FreeAlgebraBasis):
    """JT basis of the free algebra (J basis with a parameter *t*).

    Keys are ``(tuple, int)`` pairs where the tuple is a nonzero code
    and the integer tracks a power of the parameter *t*.
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a tuple or list."""
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* to a tuple key."""
        return tuple(x)

    @staticmethod
    def from_perm(perm, n):
        """Extract a JT key from *perm* if the first *n* code entries are nonzero."""
        cd = perm.code
        if len(cd) < n:
            return None
        if 0 in cd[:n]:
            return None
        return cd[:n]

    @staticmethod
    def pare_schubert(perm):
        """Extract a nonzero trimcode from *perm*, or None if it contains zeros."""
        cd = [*perm.code]
        while len(cd) > 0 and cd[-1] == 0:
            cd = cd[:-1]
        if 0 in cd:
            return None
        return tuple(cd)

    @staticmethod
    def normalize_dct(dct):
        """Normalize a word dict by collecting zeros into leading positions."""
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
        """Return a *t*-weighted ``JT``-labelled display object."""
        return sympy_Mul(JTBasis.t ** k[1], GenericPrintingTerm((k[0], "t"), "JT"))

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from JTBasis to *other_basis*."""
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

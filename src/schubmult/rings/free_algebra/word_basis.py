from functools import cache

import schubmult.rings.free_algebra as fa
from schubmult.schub_lib.permutation import uncode
from schubmult.symbolic import S, Symbol
from schubmult.utils.perm_utils import add_perm_dict

from ..schubert.schubert_ring import DSx, Sx
from ..separated_descents import SeparatedDescentsRing
from .free_algebra_basis import FreeAlgebraBasis

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)


class WordBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        return tuple(x)

    @classmethod
    def from_rc_graph(cls, rc_graph):
        return {rc_graph.length_vector(): 1}

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        return {(*key1, *key2): coeff}

    zero_monom = ()

    @classmethod
    @cache
    def coproduct(cls, key, coeff=S.One):
        if len(key) == 0:
            return {((), ()): coeff}
        if len(key) == 1:
            key = key[0]
            dct = {}
            for i in range(key + 1):
                dct[((i,), (key - i,))] = coeff
            return dct
        mid = len(key) // 2
        cp1 = cls.coproduct(key[:mid], coeff)
        cp2 = cls.coproduct(key[mid:])
        ret = {}
        for k0, v0 in cp1.items():
            for k1, v1 in cp2.items():
                ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
        return ret

    @classmethod
    @cache
    def bcoproduct(cls, key, coeff=S.One):
        if len(key) == 0:
            return {((), ()): coeff}
        if len(key) == 1:
            key = key[0]
            dct = {}
            for i in range(key + 1):
                dct[((i,) if i != 0 else (), (key - i,) if key - i != 0 else ())] = coeff
            return dct
        mid = len(key) // 2
        cp1 = cls.bcoproduct(key[:mid], coeff)
        cp2 = cls.bcoproduct(key[mid:])
        ret = {}
        for k0, v0 in cp1.items():
            for k1, v1 in cp2.items():
                ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
        return ret

    @classmethod
    def try_internal_product(cls, key1, key2, coeff=S.One):
        from sage.combinat.integer_matrices import IntegerMatrices

        bkey1 = [a + 1 for a in key1]
        bkey2 = [a + 1 for a in key2]
        mats = IntegerMatrices(bkey1, bkey2).list()

        def mat_to_key(mat):
            res = []
            for row in mat:
                for val in row:
                    if val != 0:
                        res.append(val - 1)
            return tuple(res)

        ret_dict = {}
        for mat in mats:
            key = mat_to_key(mat)
            if key not in ret_dict:
                ret_dict[key] = S.Zero
            ret_dict[key] += coeff
        return ret_dict

    @classmethod
    def internal_product(cls, key1, key2, coeff=S.One):
        if 0 in key1 or 0 in key2:
            return {}
        from sage.combinat.integer_matrices import IntegerMatrices

        mats = list(IntegerMatrices(key1, key2))

        def mat_to_key(mat):
            res = []
            for row in mat:
                for val in row:
                    if val != 0:
                        res.append(int(val))
            return tuple(res)

        ret_dict = {}
        for mat in mats:
            key = mat_to_key(mat)
            if key not in ret_dict:
                ret_dict[key] = S.Zero
            ret_dict[key] += coeff
        return ret_dict

    @classmethod
    def printing_term(cls, k):
        if all(a < 10 for a in k):
            return Symbol("[" + "".join([str(a) for a in k]) + "]")
        return Symbol("[" + " ".join([str(a) for a in k]) + "]")

    @staticmethod
    @cache
    def tup_expand(tup):
        res = splugSx([])
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            return splugSx(uncode(tup), 1)
        mid = len(tup) // 2
        return WordBasis.tup_expand(tup[:mid]) * WordBasis.tup_expand(tup[mid:])

    @staticmethod
    @cache
    def jbasis_tup_expand(tup):
        from .z_basis import ZBasis

        JB = fa.FreeAlgebra(basis=ZBasis)
        res = JB()
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            if tup[0] == 0:
                return res
            return JB(*tup)
        mid = len(tup) // 2
        return WordBasis.jbasis_tup_expand(tup[:mid]) * WordBasis.jbasis_tup_expand(tup[mid:])

    @classmethod
    def transition_schubert(cls, key):
        return dict(WordBasis.tup_expand(key))

    @classmethod
    def transition_jbasis(cls, key):
        from .j_basis import JBasis
        from .schubert_basis import SchubertBasis

        t = S.One
        JB = fa.FreeAlgebra(basis=JBasis)
        ASx = fa.FreeAlgebra(basis=SchubertBasis)
        assert all(a != 0 for a in key), "Transition from WordBasis to JBasis is only implemented for keys with no zeros."
        pangea = JB()

        for a in reversed(key):
            new_pangea = JB.from_dict({})
            for bucket, coeff in pangea.items():
                the_pieri = coeff * ASx(uncode([a])) * ASx(uncode(bucket))
                for (perm, n), v in the_pieri.items():
                    if 0 in perm.trimcode:
                        new_pangea += v * t * JB(*[a for a in perm.trimcode if a != 0])
                    else:
                        new_pangea += v * JB(*perm.trimcode)
            pangea = new_pangea

        return pangea

    @classmethod
    def transition_jtbasis(cls, key):
        from .jt_basis import JTBasis

        FA = fa.FreeAlgebra(WordBasis)
        JB = fa.FreeAlgebra(basis=JTBasis)
        res = FA.from_dict(JTBasis.normalize_dct(FA(*key)))

        ret = JB.from_dict({})

        while res != FA.zero:
            res = FA.from_dict({k: v for k, v in res.items() if v != S.Zero})
            tup = next(iter(sorted(res.keys())))
            c = res[tup]
            new_tup = tuple([int(a) for a in tup if a != 0])
            pw = len(tup) - len(new_tup)

            val = c * JB(new_tup, pw).change_basis(WordBasis)
            res -= val
            ret += c * JB(new_tup, pw)
        return dict(ret)

    @classmethod
    def transition_zbasis(cls, key):
        from .z_basis import ZBasis

        FA = fa.FreeAlgebra(WordBasis)
        JB = fa.FreeAlgebra(basis=ZBasis)
        res = FA(*key)

        ret = JB.from_dict({})
        while res != FA.zero:
            tup = next(iter(sorted(res.keys())))
            c = res[tup]
            val = c * JB(*tup).change_basis(WordBasis)
            res -= val
            ret += c * JB(*tup)
        return dict(ret)

    @classmethod
    def transition_nelementary(cls, tup):
        from sage.all import Composition

        pain = Composition(tup)
        piss = list(pain.finer())
        return {tuple([int(p) for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

    @classmethod
    def transition(cls, other_basis):
        from .j_basis import JBasis
        from .jt_basis import JTBasis
        from .nelementary_basis import NElementaryBasis
        from .schubert_basis import SchubertBasis
        from .z_basis import ZBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        if other_basis == WordBasis:
            return lambda x: {x: S.One}
        if other_basis == NElementaryBasis:
            return lambda x: cls.transition_nelementary(x)
        if other_basis == JBasis:
            return lambda x: cls.transition_jbasis(x)
        if other_basis == JTBasis:
            return lambda x: cls.transition_jtbasis(x)
        if other_basis == ZBasis:
            return lambda x: cls.transition_zbasis(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))

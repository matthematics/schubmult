from functools import cache
from itertools import zip_longest

import schubmult.rings.free_algebra as fa
from schubmult.schub_lib.perm_lib import Permutation, trimcode, uncode
from schubmult.symbolic import (
    Add,
    Integer,
    Mul,
    S,
    Symbol,
    is_of_func_type,
    sstr,
    sympify,
    sympify_sympy,
    sympy_Mul,
)
from schubmult.symmetric_polynomials import FactorialElemSym
from schubmult.utils.perm_utils import add_perm_dict, mu_A

from .abstract_schub_poly import GenericPrintingTerm
from .schubert_ring import DSx, Sx
from .separated_descents import SeparatedDescentsRing

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)


class FreeAlgebraBasis:
    @classmethod
    def is_key(cls, x): ...

    @classmethod
    def from_rc_graph(cls, rc_graph): ...

    @classmethod
    def as_key(cls, x): ...

    # @classmethod
    # def product(cls, key1, key2, coeff=S.One): ...

    # @classmethod
    # def coproduct(cls, key, coeff=S.One): ...

    @classmethod
    def transition(cls, other_basis): ...

    # @property
    # def zero_monom(cls): ...

    @classmethod
    def printing_term(cls, key): ...

    @classmethod
    def compose_transition(cls, tkeyfunc, output):
        ret = {}
        for key, v in output.items():
            ret = add_perm_dict(ret, {k: v * v0 for k, v0 in tkeyfunc(key).items()})
        return ret

    @classmethod
    def change_tensor_basis(cls, tensor_elem, basis1, basis2):
        from .tensor_ring import TensorRing

        ring1 = tensor_elem.ring.rings[0]
        ring2 = tensor_elem.ring.rings[1]
        Tring2 = TensorRing(ring1.__class__(basis=basis1), ring2.__class__(basis=basis2))
        res = Tring2.zero
        for (key1, key2), v in tensor_elem.items():
            new_elem1 = ring1(*key1).change_basis(basis1)
            new_elem2 = ring2(*key2).change_basis(basis2)
            res += v * Tring2.ext_multiply(new_elem1, new_elem2)
        return res

    # @classmethod
    # @cache
    # def coproduct(cls, key):
    #     from ._mul_utils import _tensor_product_of_dicts_first

    #     return FreeAlgebraBasis.compose_transition(
    #         lambda x: _tensor_product_of_dicts_first(SchubertBasis.transition(cls)(x[0]), SchubertBasis.transition(cls)(x[1])),
    #         FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.coproduct(y), cls.transition(SchubertBasis)(key)),
    #     )

    @classmethod
    @cache
    def coproduct(cls, key):
        from ..utils._mul_utils import _tensor_product_of_dicts_first

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(WordBasis.transition(cls)(x[0]), WordBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda y: WordBasis.coproduct(y), cls.transition(WordBasis)(key)),
        )

    @classmethod
    @cache
    def bcoproduct(cls, key):
        from ..utils._mul_utils import _tensor_product_of_dicts_first

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(WordBasis.transition(cls)(x[0]), WordBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda y: WordBasis.bcoproduct(y), cls.transition(WordBasis)(key)),
        )

    @classmethod
    @cache
    def product(cls, key1, key2, coeff=S.One):
        left = cls.transition(WordBasis)(key1)
        right = cls.transition(WordBasis)(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(WordBasis.transition(cls), WordBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret

    @classmethod
    def internal_product(cls, key1, key2, coeff=S.One):
        left = cls.transition(WordBasis)(key1)
        right = cls.transition(WordBasis)(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(WordBasis.transition(cls), WordBasis.internal_product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret


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
            # dct[((),(key,))] = coeff
            # dct[((key,), ())] = coeff
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
            # dct[((),(key,))] = coeff
            # dct[((key,), ())] = coeff
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
        # print(f"{key1=} {key2=}")
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
            # print(f"{mat=}")
            key = mat_to_key(mat)
            if key not in ret_dict:
                ret_dict[key] = S.Zero
            ret_dict[key] += coeff
        # print(f"{ret_dict=}")
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
        # return dict(WordBasis.jbasis_tup_expand(key))
        FA = fa.FreeAlgebra(WordBasis)
        JB = fa.FreeAlgebra(basis=JBasis)
        res = FA(*key)

        ret = JB.from_dict({})
        while res != FA.zero:
            tup = next(iter(sorted(res.keys())))
            c = res[tup]
            val = c * JB(*tup).change_basis(WordBasis)
            res -= val
            ret += c * JB(*tup)
        return ret

    @classmethod
    def transition_jtbasis(cls, key):
        FA = fa.FreeAlgebra(WordBasis)
        JB = fa.FreeAlgebra(basis=JTBasis)
        res = FA.from_dict(JTBasis.normalize_dct(FA(*key)))
        # t = Symbol("t")
        # res = res.kill_zero(True, t)
        # print(f"{res=}")
        # res2 = res
        # res = FA.zero

        ret = JB.from_dict({})
        # res = res2

        # print(f"{res=}")
        while res != FA.zero:
            res = FA.from_dict({k: v for k, v in res.items() if v != S.Zero})
            tup = next(iter(sorted(res.keys())))
            c = res[tup]
            new_tup = tuple([int(a) for a in tup if a != 0])
            pw = len(tup) - len(new_tup)
            # new_tup = tuple([0] * pw + new_tup)
            # dct = expand(sympify(c))
            # coeff = None
            # max_pow = 0
            # if isinstance(dct, Integer):
            #     coeff = dct
            #     val = coeff * JB(tup, 0).change_basis(WordBasis)
            #     res -= val
            #     ret += coeff * JB(tup, 0)
            # else:
            #     coeff = None
            #     # print(f"{dct=}")
            #     if isinstance(dct, Add):
            #         for term in dct.args:
            #             if isinstance(term, Mul):
            #                 if isinstance(term.args[1], Symbol) and max_pow == 0:
            #                     max_pow = 1
            #                     coeff = term.args[0] if isinstance(term.args[0], Integer) else None
            #                     if coeff is None:
            #                         raise ValueError("Coefficient is not an integer: " + sstr(term))
            #                 elif isinstance(term.args[1], Pow):
            #                     pw = int(term.args[1].args[1])
            #                     if pw > max_pow:
            #                         max_pow = pw
            #                         coeff = term.args[0] if isinstance(term.args[0], Integer) else None
            #                         if coeff is None:
            #                             raise ValueError("Coefficient is not an integer: " + sstr(term))
            #             elif isinstance(term, Symbol):
            #                 if max_pow == 0:
            #                     max_pow = 1
            #                     # coeff = t
            #                     coeff = S.One
            #             elif isinstance(term, Integer):
            #                 if max_pow == 0:
            #                     max_pow = 0
            #                     coeff = term
            #             elif isinstance(term, Pow):
            #                 pw = term.args[1]
            #                 if pw > max_pow:
            #                     max_pow = pw
            #                     #coeff = t ** pw
            #                     coeff = S.One
            #             else:
            #                 raise ValueError("Unknown term in polynomial: " + sstr(term) + " "+type(term))
            #     else:
            #         term = dct
            #         if isinstance(term, Mul):
            #             if isinstance(term.args[1], Symbol) and max_pow == 0:
            #                 max_pow = 1
            #                 coeff = term.args[0] if isinstance(term.args[0], Integer) else None
            #                 if coeff is None:
            #                     raise ValueError("Coefficient is not an integer: " + sstr(term))
            #             elif isinstance(term.args[1], Pow):
            #                 pw = int(term.args[1].args[1])
            #                 if pw > max_pow:
            #                     max_pow = pw
            #                     coeff = term.args[0] if isinstance(term.args[0], Integer) else None
            #                     if coeff is None:
            #                         raise ValueError("Coefficient is not an integer: " + sstr(term))
            #         elif isinstance(term, Symbol):
            #             if max_pow == 0:
            #                 max_pow = 1
            #                 # coeff = t
            #                 coeff = S.One
            #         elif isinstance(term, Integer):
            #             if max_pow == 0:
            #                 max_pow = 0
            #                 coeff = term
            #         elif isinstance(term, Pow):
            #             pw = term.args[1]
            #             if pw > max_pow:
            #                 max_pow = pw
            #                 #coeff = t ** pw
            #                 coeff = S.One
            #         else:
            #             raise ValueError("Unknown term in polynomial: " + sstr(term) + " "+type(term))
            #     # print(f"{coeff=}")
            val = c * JB(new_tup, pw).change_basis(WordBasis)
            res -= val
            ret += c * JB(new_tup, pw)
            # print(f"{res=}")

            # print(f"{ret=}")
        return ret

    @classmethod
    def transition_zbasis(cls, key):
        # return dict(WordBasis.jbasis_tup_expand(key))
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
        return ret

    @classmethod
    def transition_nelementary(cls, tup):
        from sage.all import Composition

        pain = Composition(tup)
        piss = list(pain.finer())
        return {tuple([int(p) for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

    @classmethod
    def transition(cls, other_basis):
        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        if other_basis == WordBasis:
            return lambda x: {x: S.One}
        if other_basis == NElementaryBasis:
            return lambda x: cls.transition_nelementary(x)
        # if other_basis == SchubertSchurBasis:
        #     return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(SchubertSchurBasis), cls.transition_schubert(x))
        if other_basis == JBasis:
            return lambda x: cls.transition_jbasis(x)
        if other_basis == JTBasis:
            return lambda x: cls.transition_jtbasis(x)
        if other_basis == ZBasis:
            return lambda x: cls.transition_zbasis(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))


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

    @staticmethod
    def pare_schubert(perm):
        cd = [*perm.code]
        while len(cd) > 0 and cd[-1] == 0:
            cd = cd[:-1]
        if 0 in cd:
            return None
        # while len(cd) > 0 and cd[0] == 0:
        #     cd = cd[1:]
        # if 0 in cd:
        #     return None
        return tuple(cd)

    # @classmethod
    # def product(cls, key1, key2, coeff=S.One):
    #     key11 = uncode([a - 1 for a in key1])
    #     key22 = uncode([a - 1 for a in key2])

    #     dct = SchubertBasis.product((key11, len(key1)), (key22, len(key2)), coeff)

    #     JB = fa.FreeAlgebra(basis=JBasis)
    #     FA = fa.FreeAlgebra()
    #     ASx = fa.FreeAlgebra(basis=SchubertBasis)
    #     #dct_out = {}
    #     dct_out = JB.zero

    #     for perm, v in dct.items():
    #         # kk = JBasis.pare_schubert(perm[0])
    #         # if 0 in perm[0].code[: perm[1]]:
    #         #     continue
    #         bento_box = ASx(*perm).change_basis(WordBasis)
    #         bento_box = ({tuple(a for a in k if a!=0): veta for k, veta in bento_box.items()})
    #         # print(f'{bento_box=}')
    #         dct_out += v * bento_box.change_basis(JBasis)
    #         # print(f"{kk=} {key1=} {key2=} {perm=}")
    #         #dct_out[tuple(kk)] = dct_out.get(tuple(kk), S.Zero) + v
    #     return dct_out

    # NDD

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "J")

    # @classmethod
    # def transition_schubert(cls, key):
    #     return dict(WordBasis.tup_expand(key))

    # @classmethod
    # def transition(cls, other_basis):
    #     # if other_basis == SchubertBasis:
    #     #     return lambda x: cls.transition_schubert(x)
    #     # if other_basis == WordBasis:
    #     #     return lambda x: {x: S.One}
    #     # if other_basis == SchubertSchurBasis:
    #     def trans(x):
    #         ASx = fa.FreeAlgebra(basis=SchubertBasis)
    #         dct = ASx(uncode(x)).change_basis(WordBasis)
    #         dct_out = {}
    #         for k, v in dct.items():
    #             k2 = tuple(a for a in k if a != 0)
    #             dct_out[k2] = dct_out.get(k2, S.Zero) + v
    #         return dct_out
    @classmethod
    def transition(cls, other_basis):
        # if other_basis == SchubertBasis:
        #     return lambda x: cls.transition_schubert(x)
        # if other_basis == WordBasis:
        #     return lambda x: {x: S.One}
        # if other_basis == SchubertSchurBasis:
        def trans(x):
            ASx = fa.FreeAlgebra(basis=SchubertBasis)
            dct = ASx(uncode(x)).change_basis(WordBasis)
            dct_out = {}
            for k, v in dct.items():
                k2 = tuple(a for a in k if a != 0)
                # if 0 in k:
                #     continue
                # k2 = k
                dct_out[k2] = dct_out.get(k2, S.Zero) + v
            return dct_out

        if other_basis == JBasis:
            return lambda x: {x: S.One}
        if other_basis == WordBasis:
            return trans

        return lambda x: FreeAlgebraBasis.compose_transition(lambda k: WordBasis.transition(other_basis)(k), trans(x))
        # return None


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
        # while len(cd) > 0 and cd[0] == 0:
        #     cd = cd[1:]
        # if 0 in cd:
        #     return None
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

    # @classmethod
    # def product(cls, key1, key2, coeff=S.One):
    #     key11 = uncode([a - 1 for a in key1])
    #     key22 = uncode([a - 1 for a in key2])

    #     dct = SchubertBasis.product((key11, len(key1)), (key22, len(key2)), coeff)

    #     JB = fa.FreeAlgebra(basis=JBasis)
    #     FA = fa.FreeAlgebra()
    #     ASx = fa.FreeAlgebra(basis=SchubertBasis)
    #     #dct_out = {}
    #     dct_out = JB.zero

    #     for perm, v in dct.items():
    #         # kk = JBasis.pare_schubert(perm[0])
    #         # if 0 in perm[0].code[: perm[1]]:
    #         #     continue
    #         bento_box = ASx(*perm).change_basis(WordBasis)
    #         bento_box = ({tuple(a for a in k if a!=0): veta for k, veta in bento_box.items()})
    #         # print(f'{bento_box=}')
    #         dct_out += v * bento_box.change_basis(JBasis)
    #         # print(f"{kk=} {key1=} {key2=} {perm=}")
    #         #dct_out[tuple(kk)] = dct_out.get(tuple(kk), S.Zero) + v
    #     return dct_out

    # NDD

    zero_monom = ()
    t = Symbol("t")

    @classmethod
    def printing_term(cls, k):
        return sympy_Mul(JTBasis.t ** k[1], GenericPrintingTerm((k[0], "t"), "JT"))

    # @classmethod
    # def transition_schubert(cls, key):
    #     return dict(WordBasis.tup_expand(key))

    # @classmethod
    # def transition(cls, other_basis):
    #     # if other_basis == SchubertBasis:
    #     #     return lambda x: cls.transition_schubert(x)
    #     # if other_basis == WordBasis:
    #     #     return lambda x: {x: S.One}
    #     # if other_basis == SchubertSchurBasis:
    #     def trans(x):
    #         ASx = fa.FreeAlgebra(basis=SchubertBasis)
    #         dct = ASx(uncode(x)).change_basis(WordBasis)
    #         dct_out = {}
    #         for k, v in dct.items():
    #             k2 = tuple(a for a in k if a != 0)
    #             dct_out[k2] = dct_out.get(k2, S.Zero) + v
    #         return dct_out
    @classmethod
    def transition(cls, other_basis):
        # if other_basis == SchubertBasis:
        #     return lambda x: cls.transition_schubert(x)
        # if other_basis == WordBasis:
        #     return lambda x: {x: S.One}
        # if other_basis == SchubertSchurBasis:
        @cache
        def trans(x2):
            # ASx = fa.FreeAlgebra(basis=SchubertBasis)
            from schubmult.rings import ASx

            if len(x2) != 2:
                raise ValueError("JTBasis transition expects a tuple of length 2, got: " + sstr(x2))
            x = x2[0]
            n = x2[1]
            # print(f"{x2=}")
            return JTBasis.normalize_dct(ASx(uncode([0] * int(n) + list(x)), n + len(x)).change_basis(WordBasis))
            # print(f"{dct=}")
            # dct_out = {}
            # for k, v in dct.items():
            #     k2 = tuple(a for a in k if a != 0)
            #     # if 0 in k:
            #     #     continue
            #     #k2 = k
            #     pw = (len(k) - len(k2))
            #     k2 = tuple([0]*pw + list(k2))
            #     assert len(k) == len(k)
            #     dct_out[k2] = dct_out.get(k2, S.Zero) + v

        if other_basis == JTBasis:
            return lambda x: {x: S.One}
        if other_basis == WordBasis:
            return trans
        return lambda x: FreeAlgebraBasis.compose_transition(lambda k: WordBasis.transition(other_basis)(k), trans(x))


class SchubertBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return (len(x) == 1 and isinstance(x[0], Permutation | list | tuple)) or (len(x) == 2 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], int))

    @classmethod
    def from_rc_graph(cls, rc_graph):
        return {(rc_graph.perm, len(rc_graph)): 1}

    @classmethod
    def as_key(cls, x):
        # print(f"{x=} {type(x)=}")
        if len(x) == 1:
            perm = Permutation(x[0])
            return (perm, 0) if len(perm.descents()) == 0 else (perm, max(perm.descents()) + 1)
        return (Permutation(x[0]), x[1])

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        return dict(coeff * splugSx(*cls.as_key(key1)) * splugSx(*cls.as_key(key2)))

    zero_monom = (Permutation([]), 0)

    @classmethod
    def skew_element(cls, w, u, n):
        from schubmult.mult.single import schubmult_py_down

        if u.inv > 0 and max(u.descents()) >= n:
            return {}
        dct = schubmult_py_down({w: S.One}, u)
        ret = {}
        for perm, v in dct.items():
            if v != S.Zero and (perm.inv == 0 or max(perm.descents()) < n):
                ret[(perm, n)] = v
        return ret

    @classmethod
    @cache
    def coproduct(cls, key):
        from ..utils._mul_utils import _tensor_product_of_dicts_first

        # print(f"{key=}")
        dct = cls.transition_word(*key)
        res = {}
        wbasis = WordBasis
        for key_word, v in dct.items():
            dct2 = wbasis.coproduct(key_word, v)
            # print(f"{dct2=}")
            for (k1, k2), v2 in dct2.items():
                dct0 = wbasis.transition_schubert(k1)
                dct1 = wbasis.transition_schubert(k2)
                res = add_perm_dict(res, {k: v0 * v2 for k, v0 in _tensor_product_of_dicts_first(dct0, dct1).items()})
        return res

    @classmethod
    def transition_schubert_schur(cls, *x):
        perm, numvars = x
        # schur(n)schubert(n)schur(m,init)schubert(n,end)
        # expansion vmu^{-1} in m schur (complement) n schubert (times w0)
        extra = len(perm) - numvars

        if extra == 0:
            return {(tuple([0] * numvars), perm, numvars): 1}
        dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        tosplit = perm * dom
        dct = Sx(tosplit).coproduct(*list(range(1, extra + 1)))
        w0 = uncode(list(range(numvars - 1, 0, -1)))
        w0s = uncode([numvars] * extra)
        dct2 = {}
        for (lambd, perm0), v in dct.items():
            perm1 = perm0 * w0
            lambd2 = tuple(trimcode(lambd * (~w0s)))
            dct2[(lambd2, perm1)] = v
        return {(k[0], *([0] * (numvars - len(k[0]))), k[1]): v for k, v in dct2.items()}

    # schubert schur => elementary multiplicative

    @classmethod
    def transition_elementary(cls, perm, numvars):
        from schubmult.utils.perm_utils import p_trans

        from .variables import genset_dict_from_expr

        mu = p_trans(list(range(numvars, 0, -1)))
        if len(mu) < len(perm) - 1:
            mu = ([numvars] * (len(perm) - 1 - len(mu))) + mu
        w0 = uncode(mu)
        dct = genset_dict_from_expr(Sx(perm * w0).as_polynomial(), Sx([]).ring.genset)
        # print(f"{dct=} ")
        # print(f"{mu=}")
        ret = {}

        for tup, v in dct.items():
            new_tup = [mu[i] - tup[i] for i in range(len(tup))]

            if len(new_tup) < len(mu):
                new_tup += mu[len(new_tup) :]
            ret[((*reversed(new_tup[-numvars + 1 :]), *sorted(new_tup[: -numvars + 1])), numvars)] = v
        return ret

    @classmethod
    def transition_separated_descents(cls, k, *x):
        perm, numvars = x
        # schur(n)schubert(n)schur(m,init)schubert(n,end)
        # expansion vmu^{-1} in m schur (complement) n schubert (times w0)
        # extra = len(perm) - numvars + k

        # if extra == 0:
        #     return {(tuple([0] * numvars), perm, numvars): 1}
        # dom = uncode([numvars] * extra + list(range(numvars - 1, 0, -1)))
        # dom = uncode(list(range(len(perm) - 1, 0, -1)))
        # tosplit = perm * dom
        # dct = Sx(tosplit).coproduct(*list(range(1, len(perm) - k + 1)))
        # w0 = uncode(list(range(k - 1, 0, -1)))
        # w0s = uncode(list(range(len(perm) - 1, k - 1, -1)))
        # dct2 = {}
        # for (perm0, perm1), v in dct.items():
        #     perm0_out = perm0 * (~w0s)
        #     perm1_out = perm1 * (w0)
        #     dct2[(perm0_out, perm1_out, numvars)] = v

        # dom0_code = perm.mul_dominant().trimcode
        dom0_code = list(range(len(perm) - 1, 0, -1))
        dom = uncode(dom0_code)
        spot = len(dom0_code) - k + 2
        tosplit = perm * (~dom)
        dct = Sx(tosplit).coproduct(*list(range(spot, len(perm))))
        w0 = uncode(mu_A(dom.code, list(range(spot - 1, len(perm)))))
        w0s = uncode(mu_A(dom.code, list(range(spot - 1))))
        # print(f"w0={w0.trimcode} w0s={w0s.trimcode} dom={dom.trimcode} spot={spot} numvars={numvars}")
        # print(f"{dct=}")
        dct2 = {}
        for (perm0, perm1), v in dct.items():
            perm0_out = perm0 * (w0)
            perm1_out = perm1 * (w0s)
            dct2[(perm1_out, perm0_out, numvars)] = v
        return dct2

    @classmethod
    def transition(cls, other_basis):
        if other_basis == SchubertBasis:
            return lambda x: {SchubertBasis.as_key(x): S.One}
        if other_basis == ElementaryBasis:
            return lambda x: cls.transition_elementary(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: cls.transition_schubert_schur(*x)
        if other_basis == WordBasis:
            return lambda x: cls.transition_word(*x)
        if other_basis.__name__ == "_SeparatedDescentsBasis":
            return lambda x: cls.transition_separated_descents(other_basis.k, *x)
        if other_basis == ElementaryBasis:
            return lambda x: cls.transition_elementary(*x)
        if other_basis == ZBasis or other_basis == JBasis or other_basis == JTBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition(WordBasis)(x))
        return None

    @classmethod
    @cache
    def transition_word(cls, perm, numvars):
        res = {}
        expr = Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis().expand()
        args = expr.args
        if not isinstance(expr, Add):
            args = [expr]
        for arg in args:
            tup = list(range(perm.inv + numvars, perm.inv, -1))
            coeff = S.One
            if is_of_func_type(sympify(arg), FactorialElemSym):
                arg = sympify_sympy(arg)
                tup[perm.inv + numvars - arg.numvars] = arg.numvars - arg.degree
            elif isinstance(arg, Mul):
                for arg0 in arg.args:
                    if is_of_func_type(sympify(arg0), FactorialElemSym):
                        arg0 = sympify_sympy(arg0)
                        tup[perm.inv + numvars - arg0.numvars] = arg0.numvars - arg0.degree
                    else:
                        coeff = Integer(arg0)
            else:
                coeff = Integer(arg)
            tup = tuple(tup)
            res[tup] = res.get(tup, S.Zero) + coeff
        return res

    @classmethod
    def printing_term(cls, k):
        return splugSx([]).ring.printing_term(k)


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
        from ..utils._mul_utils import _tensor_product_of_dicts_first

        return FreeAlgebraBasis.compose_transition(
            lambda x: _tensor_product_of_dicts_first(SchubertBasis.transition(cls)(x[0]), SchubertBasis.transition(cls)(x[1])),
            FreeAlgebraBasis.compose_transition(lambda x: SchubertBasis.coproduct(x), cls.transition_schubert(*key)),
        )
        # from ._mul_utils import _tensor_product_of_dicts_first

        # dct = cls.transition_word(*key)
        # res = {}
        # wbasis = WordBasis
        # for key_word, v in dct.items():
        #     dct2 = wbasis.coproduct(key_word, v)
        #     # print(f"{dct2=}")
        #     for (k1, k2), v2 in dct2.items():
        #         dct0 = wbasis.transition_schubert(k1)
        #         dct1 = wbasis.transition_schubert(k2)
        #         res = add_perm_dict(res, {k: v0 * v2 for k, v0 in _tensor_product_of_dicts_first(dct0, dct1).items()})
        # return res

    # def schubert_matrix(self, lambd, perm0, perm):
    # div schubert: perm * (~mu)
    # div lambd: lambd * (~mu0)
    # remaining: div(lambd * (~mu0), perm * (~mu)) on +n perm0 w9
    # multiply the +n perm0 w0

    @classmethod
    def transition_schubert(cls, lambd, perm):
        # pass
        # transition matrix, positive subsitute, Schubert times Schur
        from .schubert_ring import SingleSchubertRing
        from .variables import MaskedGeneratingSet

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
        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis), cls.transition_schubert(lambd, perm))

    @classmethod
    def transition(cls, other_basis):
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


class ElementaryBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list) and len(x) == 2 and isinstance(x[0], tuple | list) and isinstance(x[1], int)

    @classmethod
    def as_key(cls, x):
        return (x[0], x[1])

    # @classmethod
    # def product(cls, key1, key2, coeff=S.One):
    #     tup1, numvars1 = key1
    #     tup2, numvars2 = key2
    #     # so the sum of the length parititions of these is a partition
    #     length_partition1 = list(range(numvars1, 0, -1))
    #     length_partition2 = list(range(numvars2, 0, -1))
    #     if len(length_partition1) < len(tup1):
    #         length_partition1 = [numvars1] * (len(tup1) - len(length_partition1)) + length_partition1
    #     if len(length_partition2) < len(tup2):
    #         length_partition2 = [numvars2] * (len(tup2) - len(length_partition2)) + length_partition2
    #     from itertools import zip_longest
    #     total_length_partition = [a + b for a, b in zip_longest(length_partition1, length_partition2, fillvalue=0)]

    #     #return dict(coeff * splugSx(*cls.as_key(key1)) * splugSx(*cls.as_key(key2)))

    zero_monom = ((), 0)

    @classmethod
    def transition(cls, other_basis):
        if other_basis == cls:
            return lambda x: {x: S.One}
        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(*x)
        if other_basis == SchubertSchurBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_schubert_schur(*y), cls.transition_schubert(*x))
        if other_basis == WordBasis:
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_word(*y), cls.transition_schubert(*x))
        if other_basis.__name__ == "_SeparatedDescentsBasis":
            return lambda x: FreeAlgebraBasis.compose_transition(lambda y: SchubertBasis.transition_separated_descents(other_basis.k, *y), cls.transition_schubert(*x))
        return None

    @classmethod
    def transition_schubert(cls, tup, numvars):
        from schubmult.abc import x
        from schubmult.symbolic import prod

        from .poly_lib import monom_sym

        mu = list(range(numvars, 0, -1))
        if len(mu) < len(tup):
            mu = [*([numvars] * (len(tup) - len(mu))), *mu]
        tup = tuple(reversed(tup))
        flat_part = tup[: -numvars + 1]
        boink_part = tup[-numvars + 1 :]
        painted_bagel = Sx([]).ring.zero
        pickles = [mu[i] - flat_part[len(flat_part) - 1 - i] for i in range(len(flat_part))]
        painted_bagel = Sx.from_expr(monom_sym(pickles, len(flat_part), Sx([]).ring.genset))

        painted_bagel *= prod([x[i + 1] ** (mu[i] - boink_part[i - len(flat_part)]) for i in range(len(flat_part), len(mu))])
        w0 = ~uncode(mu)
        monom = {}
        for k, v in painted_bagel.items():
            if (k * w0).inv != w0.inv - k.inv:
                raise Exception
            monom[(k * w0, numvars)] = v
        return dict(monom)

    # @classmethod
    # def transition_schubert(cls, tup, numvars):
    #     from schubmult.abc import x
    #     from schubmult.symbolic import prod

    #     from .poly_lib import monom_sym

    #     mu = list(range(numvars, 0, -1))
    #     if len(mu) < len(tup):
    #         mu = [*([numvars] * (len(tup) - len(mu))), *mu]

    #     flat_part = tup[:numvars - 1]
    #     boink_part = tup[numvars - 1 :]
    #     painted_bagel = Sx([]).ring.zero
    #     pickles = [mu[i] - flat_part[len(flat_part) - 1 - i] for i in range(len(flat_part))]
    #     painted_bagel = Sx.from_expr(monom_sym(pickles, len(flat_part), Sx([]).ring.genset))

    #     painted_bagel *= prod([x[i + 1] ** (mu[i] - boink_part[len(mu) - 1 - i - len(flat_part)]) for i in range(len(flat_part), len(mu))])
    #     w0 = ~uncode(mu)
    #     monom = {}
    #     for k, v in painted_bagel.items():
    #         if (k * w0).inv != w0.inv - k.inv:
    #             raise Exception
    #         monom[(k * w0, numvars)] = v
    #     return dict(monom)

    @classmethod
    def printing_term(cls, k):
        from .abstract_schub_poly import GenericPrintingTerm

        return GenericPrintingTerm(k, "Elem")


class _SeparatedDescentsBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return len(x) == 3 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], Permutation | list | tuple) and isinstance(x[2], int) and (len(x[1]) <= x[2])

    @classmethod
    def as_key(cls, x):
        return (Permutation(x[0]), Permutation(x[1]), x[2])

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        left = cls.transition_schubert(*key1)
        right = cls.transition_schubert(*key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(SchubertBasis.transition(cls), SchubertBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret

    # @classmethod
    # @cache
    # def coproduct(cls, key):
    #     ...

    @classmethod
    def transition_schubert(cls, perm0, perm1, numvars):
        from .schubert_ring import SingleSchubertRing
        from .variables import MaskedGeneratingSet

        mu1_code = perm0.mul_dominant().trimcode
        # print(f"{mu1_code=}")
        mu2_code = list(range(cls.k - 1, 0, -1))
        # print(f"{mu2_code=}")
        mu_code = [a + b for a, b in zip_longest(mu1_code, mu2_code, fillvalue=0)]

        mu = uncode(mu_code)
        mu1 = uncode(mu1_code)
        mu2 = uncode(mu2_code)

        mul1 = perm0 * (~mu1)
        # print(f"{mul1=}")
        mul2 = perm1 * (~mu2)
        # print(f"{mul2=}")
        shifted_ring = SingleSchubertRing(MaskedGeneratingSet(Sx([]).ring.genset, list(range(1, len((~mu).trimcode) - cls.k + 2))))
        start_schub = Sx(mul1)
        # print(f"{start_schub=}")
        start_schub *= shifted_ring(mul2).in_SEM_basis()
        # print(f"{start_schub=} after mul")
        return {(k1 * mu, numvars): v for k1, v in start_schub.items()}

    @classmethod
    def transition_word(cls, perm0, perm1, n):
        return FreeAlgebraBasis.compose_transition(SchubertBasis.transition(WordBasis), cls.transition_schubert(perm0, perm1, n))

    @classmethod
    def transition(cls, other_basis):
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


# class ZNElementaryBasis(FreeAlgebraBasis):
#     @classmethod
#     def is_key(cls, x):
#         return isinstance(x, tuple | list)

#     @classmethod
#     def as_key(cls, x):
#         return tuple(x)

#     @classmethod
#     def product(cls, key1, key2, coeff=S.One):
#         return {(*key1, *key2): coeff}

#     zero_monom = ()

#     # @classmethod
#     # @cache
#     # def coproduct(cls, key, coeff=S.One):
#     #     if len(key) == 0:
#     #         return {((), ()): coeff}
#     #     if len(key) == 1:
#     #         key = key[0]
#     #         dct = {}
#     #         # dct[((),(key,))] = coeff
#     #         # dct[((key,), ())] = coeff
#     #         for i in range(key + 1):
#     #             dct[((i,), (key - i,))] = coeff
#     #         return dct
#     #     mid = len(key) // 2
#     #     cp1 = cls.coproduct(key[:mid], coeff)
#     #     cp2 = cls.coproduct(key[mid:])
#     #     ret = {}
#     #     for k0, v0 in cp1.items():
#     #         for k1, v1 in cp2.items():
#     #             ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
#     #     return ret

#     # @classmethod
#     # def try_internal_product(cls, key1, key2, coeff=S.One):
#     #     from sage.combinat.integer_matrices import IntegerMatrices

#     #     bkey1 = [a + 1 for a in key1]
#     #     bkey2 = [a + 1 for a in key2]
#     #     mats = IntegerMatrices(bkey1, bkey2).list()

#     #     def mat_to_key(mat):
#     #         res = []
#     #         for row in mat:
#     #             for val in row:
#     #                 if val != 0:
#     #                     res.append(val - 1)
#     #         return tuple(res)

#     #     ret_dict = {}
#     #     for mat in mats:
#     #         key = mat_to_key(mat)
#     #         if key not in ret_dict:
#     #             ret_dict[key] = S.Zero
#     #         ret_dict[key] += coeff
#     #     return ret_dict

#     # @classmethod
#     # def internal_product(cls, key1, key2, coeff=S.One):
#     #     if 0 in key1 or 0 in key2:
#     #         return {}
#     #     # print(f"{key1=} {key2=}")
#     #     from sage.combinat.integer_matrices import IntegerMatrices

#     #     mats = list(IntegerMatrices(key1, key2))

#     #     def mat_to_key(mat):
#     #         res = []
#     #         for row in mat:
#     #             for val in row:
#     #                 if val != 0:
#     #                     res.append(int(val))
#     #         return tuple(res)

#     #     ret_dict = {}
#     #     for mat in mats:
#     #         # print(f"{mat=}")
#     #         key = mat_to_key(mat)
#     #         if key not in ret_dict:
#     #             ret_dict[key] = S.Zero
#     #         ret_dict[key] += coeff
#     #     # print(f"{ret_dict=}")
#     #     return ret_dict

#     @classmethod
#     def printing_term(cls, k):
#         # if all(a < 10 for a in k):
#         #     return GenericPrintingTerm("L" + "".join([str(a) for a in k]) + "]")
#         # return Symbol("L[" + " ".join([str(a) for a in k]) + "]")
#         return GenericPrintingTerm(k, "LZ")

#     # @staticmethod
#     # @cache
#     # def tup_expand(tup):
#     #     res = splugSx([])
#     #     if len(tup) == 0:
#     #         return res
#     #     if len(tup) == 1:
#     #         return splugSx(uncode(tup), 1)
#     #     mid = len(tup) // 2
#     #     return WordBasis.tup_expand(tup[:mid]) * WordBasis.tup_expand(tup[mid:])

#     # @staticmethod
#     # @cache
#     # def jbasis_tup_expand(tup):
#     #     JB = fa.FreeAlgebra(basis=ZBasis)
#     #     res = JB()
#     #     if len(tup) == 0:
#     #         return res
#     #     if len(tup) == 1:
#     #         if tup[0] == 0:
#     #             return res
#     #         return JB(*tup)
#     #     mid = len(tup) // 2
#     #     return WordBasis.jbasis_tup_expand(tup[:mid]) * WordBasis.jbasis_tup_expand(tup[mid:])

#     # @classmethod
#     # def transition_schubert(cls, key):
#     #     return dict(WordBasis.tup_expand(key))

#     # @classmethod
#     # def transition_jbasis(cls, key):
#     #     # return dict(WordBasis.jbasis_tup_expand(key))
#     #     FA = fa.FreeAlgebra(WordBasis)
#     #     JB = fa.FreeAlgebra(basis=ZBasis)
#     #     res = FA(*key)

#     #     ret = JB.from_dict({})
#     #     while res != FA.zero:
#     #         tup = next(iter(sorted(res.keys())))
#     #         c = res[tup]
#     #         val = c * JB(*tup).change_basis(WordBasis)
#     #         res -= val
#     #         ret += c * JB(*tup)
#     #     return ret

#     @classmethod
#     def transition_word(cls, tup):
#         from sage.all import Composition
#         pain = Composition([a + 1 for a in tup])
#         piss = list(pain.finer())
#         return {tuple([int(p) - 1 for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

#     @classmethod
#     def transition(cls, other_basis):
#         # if other_basis == SchubertBasis:
#         #     return lambda x: cls.transition_schubert(x)
#         # if other_basis == WordBasis:
#         #     return lambda x: {x: S.One}
#         # if other_basis == SchubertSchurBasis:
#         return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(x))
#         # if other_basis == ZBasis:
#         #     return lambda x: cls.transition_jbasis(x)
#         # return None


class NElementaryBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        return tuple(x)

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        return {(*key1, *key2): coeff}

    zero_monom = ()

    # @classmethod
    # @cache
    # def coproduct(cls, key, coeff=S.One):
    #     if len(key) == 0:
    #         return {((), ()): coeff}
    #     if len(key) == 1:
    #         key = key[0]
    #         dct = {}
    #         # dct[((),(key,))] = coeff
    #         # dct[((key,), ())] = coeff
    #         for i in range(key + 1):
    #             dct[((i,), (key - i,))] = coeff
    #         return dct
    #     mid = len(key) // 2
    #     cp1 = cls.coproduct(key[:mid], coeff)
    #     cp2 = cls.coproduct(key[mid:])
    #     ret = {}
    #     for k0, v0 in cp1.items():
    #         for k1, v1 in cp2.items():
    #             ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
    #     return ret

    # @classmethod
    # def try_internal_product(cls, key1, key2, coeff=S.One):
    #     from sage.combinat.integer_matrices import IntegerMatrices

    #     bkey1 = [a + 1 for a in key1]
    #     bkey2 = [a + 1 for a in key2]
    #     mats = IntegerMatrices(bkey1, bkey2).list()

    #     def mat_to_key(mat):
    #         res = []
    #         for row in mat:
    #             for val in row:
    #                 if val != 0:
    #                     res.append(val - 1)
    #         return tuple(res)

    #     ret_dict = {}
    #     for mat in mats:
    #         key = mat_to_key(mat)
    #         if key not in ret_dict:
    #             ret_dict[key] = S.Zero
    #         ret_dict[key] += coeff
    #     return ret_dict

    # @classmethod
    # def internal_product(cls, key1, key2, coeff=S.One):
    #     if 0 in key1 or 0 in key2:
    #         return {}
    #     # print(f"{key1=} {key2=}")
    #     from sage.combinat.integer_matrices import IntegerMatrices

    #     mats = list(IntegerMatrices(key1, key2))

    #     def mat_to_key(mat):
    #         res = []
    #         for row in mat:
    #             for val in row:
    #                 if val != 0:
    #                     res.append(int(val))
    #         return tuple(res)

    #     ret_dict = {}
    #     for mat in mats:
    #         # print(f"{mat=}")
    #         key = mat_to_key(mat)
    #         if key not in ret_dict:
    #             ret_dict[key] = S.Zero
    #         ret_dict[key] += coeff
    #     # print(f"{ret_dict=}")
    #     return ret_dict

    @classmethod
    def printing_term(cls, k):
        # if all(a < 10 for a in k):
        #     return GenericPrintingTerm("L" + "".join([str(a) for a in k]) + "]")
        # return Symbol("L[" + " ".join([str(a) for a in k]) + "]")
        return GenericPrintingTerm(k, "L")

    # @staticmethod
    # @cache
    # def tup_expand(tup):
    #     res = splugSx([])
    #     if len(tup) == 0:
    #         return res
    #     if len(tup) == 1:
    #         return splugSx(uncode(tup), 1)
    #     mid = len(tup) // 2
    #     return WordBasis.tup_expand(tup[:mid]) * WordBasis.tup_expand(tup[mid:])

    # @staticmethod
    # @cache
    # def jbasis_tup_expand(tup):
    #     JB = fa.FreeAlgebra(basis=ZBasis)
    #     res = JB()
    #     if len(tup) == 0:
    #         return res
    #     if len(tup) == 1:
    #         if tup[0] == 0:
    #             return res
    #         return JB(*tup)
    #     mid = len(tup) // 2
    #     return WordBasis.jbasis_tup_expand(tup[:mid]) * WordBasis.jbasis_tup_expand(tup[mid:])

    # @classmethod
    # def transition_schubert(cls, key):
    #     return dict(WordBasis.tup_expand(key))

    # @classmethod
    # def transition_jbasis(cls, key):
    #     # return dict(WordBasis.jbasis_tup_expand(key))
    #     FA = fa.FreeAlgebra(WordBasis)
    #     JB = fa.FreeAlgebra(basis=ZBasis)
    #     res = FA(*key)

    #     ret = JB.from_dict({})
    #     while res != FA.zero:
    #         tup = next(iter(sorted(res.keys())))
    #         c = res[tup]
    #         val = c * JB(*tup).change_basis(WordBasis)
    #         res -= val
    #         ret += c * JB(*tup)
    #     return ret

    @classmethod
    def transition_word(cls, tup):
        from sage.all import Composition

        pain = Composition(tup)
        piss = list(pain.finer())
        return {tuple([int(p) for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

    @classmethod
    def transition(cls, other_basis):
        # if other_basis == SchubertBasis:
        #     return lambda x: cls.transition_schubert(x)
        # if other_basis == WordBasis:
        #     return lambda x: {x: S.One}
        # if other_basis == SchubertSchurBasis:
        return lambda x: FreeAlgebraBasis.compose_transition(WordBasis.transition(other_basis), cls.transition_word(x))


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
        # while len(cd) > 0 and cd[0] == 0:
        #     cd = cd[1:]
        # if 0 in cd:
        #     return None
        return tuple(cd)

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        key11 = uncode([a - 1 for a in key1])
        key22 = uncode([a - 1 for a in key2])

        dct = SchubertBasis.product((key11, len(key1)), (key22, len(key2)), coeff)

        JB = fa.FreeAlgebra(basis=ZBasis)
        # FA = fa.FreeAlgebra()
        # ASx = fa.FreeAlgebra(basis=SchubertBasis)
        # dct_out = {}
        dct_out = JB.zero

        for (perm, n), v in dct.items():
            # kk = ZBasis.pare_schubert(perm[0])
            # if 0 in perm[0].code[: perm[1]]:
            #     continue
            cd = [*perm.code]
            tup = [1] * n
            for i, a in enumerate(cd[:n]):
                tup[i] = a + 1
            dct_out[tuple(tup)] = v
            # print(f"{kk=} {key1=} {key2=} {perm=}")
            # dct_out[tuple(kk)] = dct_out.get(tuple(kk), S.Zero) + v
        return dct_out

    # NDD

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "Z")

    # @classmethod
    # def transition_schubert(cls, key):
    #     return dict(WordBasis.tup_expand(key))

    # @classmethod
    # def transition(cls, other_basis):
    #     # if other_basis == SchubertBasis:
    #     #     return lambda x: cls.transition_schubert(x)
    #     # if other_basis == WordBasis:
    #     #     return lambda x: {x: S.One}
    #     # if other_basis == SchubertSchurBasis:
    #     def trans(x):
    #         ASx = fa.FreeAlgebra(basis=SchubertBasis)
    #         dct = ASx(uncode(x)).change_basis(WordBasis)
    #         dct_out = {}
    #         for k, v in dct.items():
    #             k2 = tuple(a for a in k if a != 0)
    #             dct_out[k2] = dct_out.get(k2, S.Zero) + v
    #         return dct_out
    @classmethod
    def transition(cls, other_basis):
        def trans(x):
            from schubmult.rings import ASx

            dct = ASx(uncode(x)).change_basis(WordBasis)
            dct_out = {}
            for k, v in dct.items():
                # k2 = tuple(a for a in k if a != 0)
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
        # return None


# class ZBasis(FreeAlgebraBasis):
#     @classmethod
#     def is_key(cls, x):
#         return isinstance(x, tuple | list)

#     @classmethod
#     def as_key(cls, x):
#         return tuple(x)

#     @staticmethod
#     def from_perm(perm, n):
#         cd = perm.code
#         if len(cd) < n:
#             return None
#         if 0 in cd[:n]:
#             return None
#         return cd[:n]

#     @staticmethod
#     def pare_schubert(perm):
#         cd = [*perm.code]
#         while len(cd) > 0 and cd[-1] == 0:
#             cd = cd[:-1]
#         if 0 in cd:
#             return None
#         # while len(cd) > 0 and cd[0] == 0:
#         #     cd = cd[1:]
#         # if 0 in cd:
#         #     return None
#         return tuple(cd)

#     # @classmethod
#     # def product(cls, key1, key2, coeff=S.One):
#     #     key11 = uncode([a - 1 for a in key1])
#     #     key22 = uncode([a - 1 for a in key2])

#     #     dct = SchubertBasis.product((key11, len(key1)), (key22, len(key2)), coeff)

#     #     JB = fa.FreeAlgebra(basis=ZBasis)
#     #     FA = fa.FreeAlgebra()
#     #     ASx = fa.FreeAlgebra(basis=SchubertBasis)
#     #     #dct_out = {}
#     #     dct_out = JB.zero

#     #     for perm, v in dct.items():
#     #         # kk = ZBasis.pare_schubert(perm[0])
#     #         # if 0 in perm[0].code[: perm[1]]:
#     #         #     continue
#     #         bento_box = ASx(*perm).change_basis(WordBasis)
#     #         bento_box = FA.from_dict({tuple(a for a in k if a!=0): veta for k, veta in bento_box.items()})
#     #         # print(f'{bento_box=}')
#     #         dct_out += v * bento_box.change_basis(ZBasis)
#     #         # print(f"{kk=} {key1=} {key2=} {perm=}")
#     #         #dct_out[tuple(kk)] = dct_out.get(tuple(kk), S.Zero) + v
#     #     return dct_out

#     # NDD

#     zero_monom = ()

#     @classmethod
#     def printing_term(cls, k):
#         return GenericPrintingTerm(k, "Z")

#     # @classmethod
#     # def transition_schubert(cls, key):
#     #     return dict(WordBasis.tup_expand(key))

#     # @classmethod
#     # def transition(cls, other_basis):
#     #     # if other_basis == SchubertBasis:
#     #     #     return lambda x: cls.transition_schubert(x)
#     #     # if other_basis == WordBasis:
#     #     #     return lambda x: {x: S.One}
#     #     # if other_basis == SchubertSchurBasis:
#     #     def trans(x):
#     #         ASx = fa.FreeAlgebra(basis=SchubertBasis)
#     #         dct = ASx(uncode(x)).change_basis(WordBasis)
#     #         dct_out = {}
#     #         for k, v in dct.items():
#     #             k2 = tuple(a for a in k if a != 0)
#     #             dct_out[k2] = dct_out.get(k2, S.Zero) + v
#     #         return dct_out
#     @classmethod
#     def transition(cls, other_basis):
#         # if other_basis == SchubertBasis:
#         #     return lambda x: cls.transition_schubert(x)
#         # if other_basis == WordBasis:
#         #     return lambda x: {x: S.One}
#         # if other_basis == SchubertSchurBasis:
#         def trans(x):
#             ASx = fa.FreeAlgebra(basis=SchubertBasis)
#             dct = ASx(uncode(x)).change_basis(WordBasis)
#             dct_out = {}
#             for k, v in dct.items():
#                 #k2 = tuple(a for a in k if a != 0)
#                 if 0 in k:
#                     continue
#                 k2 = k
#                 dct_out[k2] = dct_out.get(k2, S.Zero) + v
#             return dct_out

#         if other_basis == ZBasis:
#             return lambda x: {x: S.One}
#         if other_basis == WordBasis:
#             return trans

#         return lambda x: FreeAlgebraBasis.compose_transition(lambda k: WordBasis.transition(other_basis)(k), trans(x))
#         # return None

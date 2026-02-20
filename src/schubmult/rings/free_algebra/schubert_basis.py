from functools import cache

from schubmult.schub_lib.permutation import Permutation, uncode
from schubmult.symbolic import Add, Integer, Mul, S, is_of_func_type, sympify, sympify_sympy
from schubmult.symmetric_polynomials import FactorialElemSym
from schubmult.utils.perm_utils import add_perm_dict, mu_A

from ..schubert.schubert_ring import DSx, Sx
from ..schubert.separated_descents import SeparatedDescentsRing
from .free_algebra_basis import FreeAlgebraBasis

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)


class SchubertBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return (len(x) == 1 and isinstance(x[0], Permutation | list | tuple)) or (len(x) == 2 and isinstance(x[0], Permutation | list | tuple) and isinstance(x[1], int))

    @classmethod
    def from_rc_graph(cls, rc_graph):
        return {(rc_graph.perm, len(rc_graph)): 1}

    @classmethod
    def as_key(cls, x):
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
        from ...utils._mul_utils import _tensor_product_of_dicts_first
        from .word_basis import WordBasis

        dct = cls.transition_word(*key)
        res = {}
        wbasis = WordBasis
        for key_word, v in dct.items():
            dct2 = wbasis.coproduct(key_word, v)
            for (k1, k2), v2 in dct2.items():
                dct0 = wbasis.transition_schubert(k1)
                dct1 = wbasis.transition_schubert(k2)
                res = add_perm_dict(res, {k: v0 * v2 for k, v0 in _tensor_product_of_dicts_first(dct0, dct1).items()})
        return res

    @classmethod
    def transition_schubert_schur(cls, *x):
        perm, numvars = x
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
            lambd2 = tuple((lambd * (~w0s)).trimcode)
            dct2[(lambd2, perm1)] = v
        return {(k[0], *([0] * (numvars - len(k[0]))), k[1]): v for k, v in dct2.items()}

    @classmethod
    def transition_elementary(cls, perm, numvars):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr
        from schubmult.utils.perm_utils import p_trans

        mu = p_trans(list(range(numvars, 0, -1)))
        if len(mu) < len(perm) - 1:
            mu = ([numvars] * (len(perm) - 1 - len(mu))) + mu
        w0 = uncode(mu)
        dct = genset_dict_from_expr(Sx(perm * w0).as_polynomial(), Sx([]).ring.genset)
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
        dom0_code = list(range(len(perm) - 1, 0, -1))
        dom = uncode(dom0_code)
        spot = len(dom0_code) - k + 2
        tosplit = perm * (~dom)
        dct = Sx(tosplit).coproduct(*list(range(spot, len(perm))))
        w0 = uncode(mu_A(dom.code, list(range(spot - 1, len(perm)))))
        w0s = uncode(mu_A(dom.code, list(range(spot - 1))))
        dct2 = {}
        for (perm0, perm1), v in dct.items():
            perm0_out = perm0 * (w0)
            perm1_out = perm1 * (w0s)
            dct2[(perm1_out, perm0_out, numvars)] = v
        return dct2

    @classmethod
    def transition_jbasis(cls, perm, n):
        from .j_basis import JBasis
        from .word_basis import WordBasis

        t = S.One
        if len(perm.trimcode) < n:
            return FreeAlgebraBasis.compose_transition(WordBasis.transition(JBasis), cls.transition_word(perm, n))
        if 0 not in perm.trimcode:
            return {tuple(perm.trimcode): S.One}
        leading_zeros = 0
        codecode = [*perm.trimcode]
        for a in perm.trimcode:
            if a == 0:
                leading_zeros += 1
                codecode = codecode[1:]
            else:
                break
        if 0 not in codecode:
            return {tuple(codecode): t ** leading_zeros}
        if leading_zeros > 0:
            return {k: v * t**leading_zeros for k, v in FreeAlgebraBasis.compose_transition(WordBasis.transition_jbasis, cls.transition_word(uncode(codecode), n - leading_zeros)).items()}
        return FreeAlgebraBasis.compose_transition(WordBasis.transition_jbasis, cls.transition_word(perm, n))

    @classmethod
    def transition(cls, other_basis):
        from .elementary_basis import ElementaryBasis
        from .j_basis import JBasis
        from .jt_basis import JTBasis
        from .schubert_schur_basis import SchubertSchurBasis
        from .word_basis import WordBasis
        from .z_basis import ZBasis

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
        if other_basis == ZBasis or other_basis == JTBasis or other_basis == JBasis:
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

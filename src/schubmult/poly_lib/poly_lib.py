from functools import cache, cached_property

import symengine
from symengine import Mul, Pow, sympify

from .variables import GeneratingSet, is_indexed


#Indexed._sympystr = lambda x, p: f"{p.doprint(x.args[0])}_{x.args[1]}"
def expand(val):
    return symengine.expand(val)


class _gvars:
    @cached_property
    def n(self):
        return 100

    # @cached_property
    # def fvar(self):
    #     return 100

    @cached_property
    def var1(self):
        return GeneratingSet("x")

    @cached_property
    def var2(self):
        return GeneratingSet("y")

    @cached_property
    def var3(self):
        return GeneratingSet("z")

    @cached_property
    def var_r(self):
        return GeneratingSet("r")

    @cached_property
    def var_g1(self):
        return GeneratingSet("y")

    @cached_property
    def var_g2(self):
        return GeneratingSet("z")

    @cached_property
    def q_var(self):
        return GeneratingSet("q")


zero = sympify(0)

_vars = _gvars()


def elem_sym_func(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2):
    newk = k - udiff
    if newk < vdiff:
        return zero
    if newk == vdiff:
        return one
    yvars = []
    for j in range(min(len(u1), k)):
        if u1[j] == u2[j]:
            yvars += [varl1[u2[j]]]
    for j in range(len(u1), min(k, len(u2))):
        if u2[j] == j + 1:
            yvars += [varl1[u2[j]]]
    for j in range(len(u2), k):
        yvars += [varl1[j + 1]]
    zvars = [varl2[i] for i in call_zvars(v1, v2, k, i)]
    return elem_sym_poly(newk - vdiff, newk, yvars, zvars)


def elem_sym_func_q(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2):
    newk = k - udiff
    if newk < vdiff:
        return zero
    if newk == vdiff:
        return one
    yvars = []
    mlen = max(len(u1), len(u2))
    u1 = [*u1] + [a + 1 for a in range(len(u1), mlen)]
    u2 = [*u2] + [a + 1 for a in range(len(u2), mlen)]
    for j in range(min(len(u1), k)):
        if u1[j] == u2[j]:
            yvars += [varl1[u2[j]]]
    for j in range(len(u1), min(k, len(u2))):
        if u2[j] == j + 1:
            yvars += [varl1[u2[j]]]
    for j in range(len(u2), k):
        yvars += [varl1[j + 1]]
    zvars = [varl2[a] for a in call_zvars(v1, v2, k, i)]
    return elem_sym_poly(newk - vdiff, newk, yvars, zvars)


one = sympify(1)


def elem_sym_poly_q(p, k, varl1, varl2, q_var=_vars.q_var):
    if p == 0 and k >= 0:
        return one
    if p < 0 or p > k:
        return zero
    return (
        (varl1[k - 1] - varl2[k - p]) * elem_sym_poly_q(p - 1, k - 1, varl1, varl2, q_var)
        + elem_sym_poly_q(p, k - 1, varl1, varl2, q_var)
        + q_var[k - 1] * elem_sym_poly_q(p - 2, k - 2, varl1, varl2, q_var)
    )


def elem_sym_poly(p, k, varl1, varl2, xstart=0, ystart=0):
    if p > k:
        return zero
    if p == 0:
        return one
    if p == 1:
        res = varl1[xstart] - varl2[ystart]
        for i in range(1, k):
            res += varl1[xstart + i] - varl2[ystart + i]
        return res
    if p == k:
        res = (varl1[xstart] - varl2[ystart]) * (varl1[xstart + 1] - varl2[ystart])
        for i in range(2, k):
            res *= varl1[i + xstart] - varl2[ystart]
        return res
    mid = k // 2
    xsm = xstart + mid
    ysm = ystart + mid
    kmm = k - mid
    res = elem_sym_poly(p, mid, varl1, varl2, xstart, ystart) + elem_sym_poly(
        p,
        kmm,
        varl1,
        varl2,
        xsm,
        ysm,
    )
    for p2 in range(max(1, p - kmm), min(p, mid + 1)):
        res += elem_sym_poly(p2, mid, varl1, varl2, xstart, ystart) * elem_sym_poly(
            p - p2,
            kmm,
            varl1,
            varl2,
            xsm,
            ysm - p2,
        )
    return res


@cache
def call_zvars(v1, v2, k, i):  # noqa: ARG001
    v3 = [*v2, *list(range(len(v2) + 1, i + 1))]
    return [v3[i - 1]] + [v3[j] for j in range(len(v1), len(v3)) if v3[j] != j + 1 and j != i - 1] + [v3[j] for j in range(len(v1)) if v1[j] != v3[j] and j != i - 1]


def efficient_subs(expr, subs_dict):
    subs_dict_new = {}
    expr = sympify(expr)
    for s in expr.free_symbols:
        if s in subs_dict:
            subs_dict_new[s] = subs_dict[s]
    return expr.subs(subs_dict_new)

def q_vector(q_exp, q_var=_vars.q_var):
    # qvar_list = q_var.tolist()
    ret = []

    if q_exp == 1:
        return ret
    if is_indexed(q_exp) and q_exp.base == q_var:
        i = q_exp.index
        return [0 for j in range(i - 1)] + [1]
    if isinstance(q_exp, Pow):
        qv = q_exp.args[0]
        expon = int(q_exp.args[1])
        i = qv.args[1]
        return [0 for j in range(i - 1)] + [expon]
    if isinstance(q_exp, Mul):
        for a in q_exp.args:
            v1 = q_vector(a)
            v1 += [0 for i in range(len(v1), len(ret))]
            ret += [0 for i in range(len(ret), len(v1))]
            ret = [ret[i] + v1[i] for i in range(len(ret))]
        return ret

    return None


def xreplace_genvars(poly, vars1, vars2):
    subs_dict = {}
    for s in sympify(poly).free_symbols:
        if is_indexed(s) and s.base == vars1:
            subs_dict[s] = vars2[s.index]
    return sympify(poly).xreplace(subs_dict)
    # print(f"{poly2=} {poly2.free_symbols=}")

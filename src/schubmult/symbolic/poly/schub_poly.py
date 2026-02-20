from functools import cache, cached_property

import schubmult.combinatorics.permutation as pl
import schubmult.utils.schub_lib as schub_lib
from schubmult.symbolic import Add, Mul, Pow, S, Symbol, prod, sympify

from . import variables as vv


class _gvars:
    @cached_property
    def n(self):
        return 100

    @cached_property
    def var1(self):
        return vv.GeneratingSet("x")

    @cached_property
    def var2(self):
        return vv.GeneratingSet("y")

    @cached_property
    def var3(self):
        return vv.GeneratingSet("z")

    @cached_property
    def var_r(self):
        return vv.GeneratingSet("r")

    @cached_property
    def var_g1(self):
        return vv.GeneratingSet("ॅ")

    @cached_property
    def var_g2(self):
        return vv.GeneratingSet("१")

    @cached_property
    def q_var(self):
        return vv.GeneratingSet("q")


zero = 0
one = 1

_vars = _gvars()


def sv_posify(val, var2):
    var_r = vv.GeneratingSet("r")
    subs_dict = {}
    for i in range(1, 100):
        sm = sympify(var2[1])
        for j in range(1, i):
            sm += var_r[j]
        subs_dict[sympify(var2[i])] = sm
    val = sympify(efficient_subs(sympify(val), subs_dict).simplify())
    bingle_dict = {}
    for i in range(1, len(var_r) - 1):
        bingle_dict[var_r[i]] = var2[i + 1] - var2[i]
    return val.xreplace(bingle_dict)


def act(w, poly, genset):
    if not isinstance(w, pl.Permutation):
        w = pl.Permutation(w)
    subs_dict = {}
    if not isinstance(genset, vv.GeneratingSet_base):
        genset = vv.CustomGeneratingSet(genset)
    for s in poly.free_symbols:
        if genset.index(s) != -1:
            subs_dict[s] = genset[w(genset.index(s))]
    return efficient_subs(poly, subs_dict)


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
    for j in range(k):
        if u1[j] == u2[j]:
            yvars += [varl1[u2[j]]]
    zvars = [varl2[a] for a in call_zvars(v1, v2, k, i)]
    return elem_sym_poly(newk - vdiff, newk, yvars, zvars)


def elem_sym_poly_q(p, k, varl1, varl2, q_var=_vars.q_var):
    if p == 0 and k >= 0:
        return S.One
    if p < 0 or p > k:
        return S.Zero
    return (
        (varl1[k - 1] - varl2[k - p]) * elem_sym_poly_q(p - 1, k - 1, varl1, varl2, q_var)
        + elem_sym_poly_q(p, k - 1, varl1, varl2, q_var)
        + q_var[k - 1] * elem_sym_poly_q(p - 2, k - 2, varl1, varl2, q_var)
    )


def complete_sym_poly(p, k, vrs, vrs2):
    if p == 0 and k >= 0:
        return S.One
    if p != 0 and k == 0:
        return S.Zero
    if k < 0:
        return S.Zero
    if k == 1:
        return prod([vrs[0] - vrs2[i] for i in range(p)])
    sm = 0
    mid = k // 2
    for i in range(p + 1):
        sm += complete_sym_poly(i, mid, vrs[:mid], vrs2[: mid + i - 1]) * complete_sym_poly(p - i, k - mid, vrs[mid:], vrs2[mid + i :])
    return sm


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
    return [v2[i - 1]] + [v2[j] for j in range(len(v1), len(v2) + max(0, i - len(v2))) if v2[j] != j + 1 and j != i - 1] + [v2[j] for j in range(len(v1)) if v1[j] != v2[j] and j != i - 1]


def efficient_subs(expr, subs_dict):
    subs_dict_new = {}
    expr = sympify(expr)
    for s in expr.free_symbols:
        if s in subs_dict:
            subs_dict_new[s] = subs_dict[s]
    return expr.subs(subs_dict_new)


def q_vector(q_exp, q_var=_vars.q_var):
    ret = []

    if q_exp == 1:
        return ret
    if q_var.index(q_exp) != -1:
        i = q_var.index(q_exp)
        return [0 for j in range(i - 1)] + [1]
    if isinstance(q_exp, Pow):
        qv = q_exp.args[0]
        expon = int(q_exp.args[1])
        i = q_var.index(qv)
        if i == -1:
            raise IndexError
        return [0 for j in range(i - 1)] + [expon]
    if isinstance(q_exp, Mul):
        for a in q_exp.args:
            v1 = q_vector(a)
            v1 += [0 for i in range(len(v1), len(ret))]
            ret += [0 for i in range(len(ret), len(v1))]
            ret = [ret[i] + v1[i] for i in range(len(ret))]
        return ret

    return None


def monom_sym(partition, numvars, genset):
    if numvars == 0:
        return S.One
    if numvars < 0:
        return S.Zero
    if len(partition) < numvars:
        partition = [*partition, *([0] * (numvars - len(partition)))]
    pm1 = -1
    res = S.Zero
    for i, p in enumerate(partition):
        if pm1 != p:
            pm1 = p
            res += (genset[numvars] ** p) * monom_sym(partition[:i] + partition[i + 1 :], numvars - 1, genset)
    return res


def xreplace_genvars(poly, vars1, vars2):
    subs_dict = {}
    for s in sympify(poly).free_symbols:
        if _vars.var_g1.index(s) != -1:
            subs_dict[s] = sympify(vars1[_vars.var_g1.index(s)])
        elif _vars.var_g2.index(s) != -1:
            subs_dict[s] = sympify(vars2[_vars.var_g2.index(s)])
    return sympify(poly).xreplace(subs_dict)


def divide_out_diff(poly, v1, v2):
    if hasattr(poly, "_eval_divide_out_diff"):
        return poly._eval_divide_out_diff(v1, v2)
    Mul_local = Mul
    Add_local = Add
    Pow_local = Pow
    sympify_local = sympify
    poly = sympify(poly)
    v1 = sympify(v1)
    v2 = sympify(v2)

    poly2 = poly.xreplace({v1: v2})
    if poly == poly2:
        return S.Zero
    if poly == v2:
        return S.NegativeOne
    if poly2 == v2:
        return S.One
    if isinstance(poly, Add_local):
        return sympify_local(sympify(Add_local(*[divide_out_diff(a, v1, v2) for a in poly.args])))
    if isinstance(poly, Pow_local):
        b = poly.args[0]
        dd = divide_out_diff(poly.args[0], v1, v2)
        a = poly.args[0].xreplace({v1: v2})
        return Add_local(*[Mul_local(dd, Pow_local(b, i), Pow_local(a, int(poly.args[1]) - 1 - i)) for i in range(int(poly.args[1]))])
    if isinstance(poly, Mul_local):
        current_args = [*poly.args]
        args_ret = []
        for i, arg in enumerate(poly.args):
            res = divide_out_diff(arg, v1, v2)
            if res == S.Zero:
                continue
            if res == S.One:
                args_ret += [Mul_local(*[*current_args[:i], *current_args[i + 1 :]])]
            else:
                args_ret += [Mul_local(*[*current_args[:i], res, *current_args[i + 1 :]])]
            current_args[i] = current_args[i].xreplace({v1: v2})
        return sympify_local(sympify(Add_local(*args_ret)))
    raise ValueError(f"Expected Expr but got {type(poly)}")


def split_up(poly, v1, v2):
    try:
        return (sympify(poly).xreplace({sympify(v1): sympify(v2)}), (sympify(v1 - v2), divide_out_diff(poly, v1, v2)))
    except Exception:
        return (sympify(poly).xreplace({sympify(v1): sympify(v2)}), (sympify(v1 - v2), divide_out_diff(poly, v1, v2)))


def perm_act(val, i, var2=None):
    subsdict = {var2[i]: var2[i + 1], var2[i + 1]: var2[i]}
    return sympify(val).subs(subsdict)


def elem_func_func(k, i, v1, v2, vdiff, varl1, varl2, elem_func):
    newk = k
    if newk < vdiff:
        return 0
    if newk == vdiff:
        return 1
    zvars = [varl2[i] for i in call_zvars(v1, v2, k, i)]
    return elem_func(newk - vdiff, newk, varl1[1:], zvars)


def elem_func_func_mul(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2, elem_func):
    newk = k - udiff
    if newk < vdiff:
        return 0
    if newk == vdiff:
        return 1
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
    return elem_func(newk - vdiff, newk, yvars, zvars)


def schubpoly_from_elems(v, var_x=None, var_y=None, elem_func=None, mumu=None):
    if mumu:
        th = mumu.code
        mu = mumu
    else:
        th = (~pl.Permutation(v)).strict_theta()
        mu = pl.uncode(th)
    vmu = pl.Permutation(v) * mu
    if len(th) == 0:
        return elem_func(0, 0, var_x, var_y)
    while len(th) > 0 and th[-1] == 0:
        th.pop()
    vpathdicts = schub_lib.compute_vpathdicts(th, vmu)
    vpathsums = {pl.Permutation([1, 2]): elem_func(0, 0, var_x, var_y)}
    for index in range(len(th)):
        mx_th = 0
        newpathsums = {}
        for vp in vpathdicts[index]:
            for v2, vdiff, s in vpathdicts[index][vp]:
                mx_th = max(mx_th, th[index] - vdiff)
        for v in vpathdicts[index]:
            sumval = vpathsums.get(v, S.Zero)
            if sumval == S.Zero:
                continue
            for v2, vdiff, s in vpathdicts[index][v]:
                newpathsums[v2] = newpathsums.get(
                    v2,
                    S.Zero,
                ) + s * sumval * elem_func_func(
                    th[index],
                    index + 1,
                    v,
                    v2,
                    vdiff,
                    var_x,
                    var_y,
                    elem_func=elem_func,
                )
        vpathsums = newpathsums
    return vpathsums.get(vmu, S.Zero)


def schubpoly_classical_from_elems(v, var_x=None, var_y=None, elem_func=None):
    th = (~pl.Permutation(v)).theta()
    mu = pl.uncode(th)
    vmu = pl.Permutation(v) * mu
    if len(th) == 0:
        return elem_func(0, 0, var_x, var_y)
    while len(th) > 0 and th[-1] == 0:
        th.pop()
    vpathdicts = schub_lib.compute_vpathdicts(th, vmu)
    vpathsums = {pl.Permutation([1, 2]): elem_func(0, 0, var_x, var_y)}
    for index in range(len(th)):
        mx_th = 0
        newpathsums = {}
        for vp in vpathdicts[index]:
            for v2, vdiff, s in vpathdicts[index][vp]:
                mx_th = max(mx_th, th[index] - vdiff)
        for v in vpathdicts[index]:
            sumval = vpathsums.get(v, S.Zero)
            if sumval == S.Zero:
                continue
            for v2, vdiff, s in vpathdicts[index][v]:
                newpathsums[v2] = newpathsums.get(
                    v2,
                    S.Zero,
                ) + s * sumval * elem_func_func(
                    th[index],
                    index + 1,
                    v,
                    v2,
                    vdiff,
                    var_x,
                    var_y,
                    elem_func=elem_func,
                )
        vpathsums = newpathsums
    return vpathsums.get(vmu, S.Zero)


def schubpoly(v, var2=None, var3=None, start_var=1):
    n = 0
    for j in range(len(v) - 2, -1, -1):
        if v[j] > v[j + 1]:
            n = j + 1
            break
    if n == 0:
        return 1
    lst = schub_lib.pull_out_var(n, v)
    ret = 0
    for pw, vp in lst:
        tomul = 1
        for p in pw:
            tomul *= var2[start_var + n - 1] - var3[p]
        ret += tomul * schubpoly(vp, var2, var3, start_var)
    return ret


_s = Symbol("_s")


def div_diff(poly, v1, v2):
    if hasattr(poly, "_eval_div_diff"):
        return poly._eval_div_diff(v1, v2)
    poly = sympify(poly)
    v1 = sympify(v1)
    v2 = sympify(v2)

    if v1 not in poly.free_symbols and v2 not in poly.free_symbols:
        return S.Zero
    if poly == v1:
        return S.One
    if poly == v2:
        return S.NegativeOne
    if isinstance(poly, Add):
        return Add(*[div_diff(a, v1, v2) for a in poly.args])
    if isinstance(poly, Pow):
        a = poly.args[0]
        dd = div_diff(poly.args[0], v1, v2)
        b = poly.args[0].xreplace({v1: _s}).xreplace({v2: v1}).xreplace({_s: v2})
        return Add(*[Mul(dd, Pow(b, i), Pow(a, int(poly.args[1]) - 1 - i)) for i in range(int(poly.args[1]))])
    if isinstance(poly, Mul):
        current_args = [*poly.args]
        args_ret = []
        for i, arg in enumerate(poly.args):
            res = div_diff(arg, v1, v2)
            if res == S.Zero:
                continue
            if res == S.One:
                args_ret += [Mul(*[*current_args[:i], *current_args[i + 1 :]])]
            else:
                args_ret += [Mul(*[*current_args[:i], res, *current_args[i + 1 :]])]
            current_args[i] = current_args[i].xreplace({v1: _s}).xreplace({v2: v1}).xreplace({_s: v2})
        return Add(*args_ret)
    raise ValueError(f"Expected Expr but got {type(poly)}")

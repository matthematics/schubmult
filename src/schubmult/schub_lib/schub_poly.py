import schubmult.schub_lib.perm_lib as pl
import schubmult.utils.schub_lib as schub_lib
from schubmult.rings.poly_lib import call_zvars
from schubmult.symbolic import Add, Mul, Pow, S, Symbol, sympify


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
        th = pl.code(mumu)
        mu = mumu
    else:
        th = pl.strict_theta(~pl.Permutation(v))
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
    # print(f"{v=} {var_x=} {var_y=}")
    th = pl.theta(~pl.Permutation(v))
    mu = pl.uncode(th)
    vmu = pl.Permutation(v) * mu  # permtrim(mulperm([*v], mu))
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
    # print("prantix")
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


# def skew_div_diff(u, w, poly):
#     d = -1
#     for i in range(len(w) - 1):
#         if w[i] > w[i + 1]:
#             d = i
#             break
#     d2 = -1
#     for i in range(len(u) - 1):
#         if u[i] > u[i + 1]:
#             d2 = i
#             break
#     if d == -1:
#         if d2 == -1:
#             return poly
#         return 0
#     w2 = w.swap(d, d + 1)
#     if d < len(u) - 1 and u[d] > u[d + 1]:
#         u2 = u.swap(d, d + 1)
#         return skew_div_diff(u2, w2, perm_act(poly, d + 1))
#     return skew_div_diff(u, w2, div_diff(d + 1, poly))

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
        b = poly.args[0].xreplace({v1: _s}).xreplace({v2:v1}).xreplace({_s:v2})
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
            current_args[i] = current_args[i].xreplace({v1: _s}).xreplace({v2:v1}).xreplace({_s:v2})
        return Add(*args_ret)
    raise ValueError(f"Expected Expr but got {type(poly)}")

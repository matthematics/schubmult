import sympy

import schubmult.perm_lib as pl
import schubmult.schub_lib.schub_lib as schub_lib
from schubmult.poly_lib.poly_lib import call_zvars
from schubmult.poly_lib.variables import GeneratingSet


def perm_act(val, i, var2=None):
    subsdict = {var2[i]: var2[i + 1], var2[i + 1]: var2[i]}
    return sympy.sympify(val).subs(subsdict)


def div_diff(i, poly, var2=None):
    return sympy.sympify(
        sympy.div(sympy.sympify(poly - perm_act(poly, i)), sympy.sympify(var2[i] - var2[i + 1]))[0],
    )


def elem_func_func(k, i, v1, v2, vdiff, varl1, varl2, elem_func):
    newk = k
    if newk < vdiff:
        return 0
    if newk == vdiff:
        return 1
    zvars = [varl2[i] for i in call_zvars(v1, v2, k, i)]
    return elem_func(newk - vdiff, newk, varl1[1:], zvars)


def schubpoly_from_elems(v, var_x=None, var_y=None, elem_func=None, mumu=None):
    if mumu:
        # print(pl.code(mumu))
        th = pl.code(mumu)
        # print(f"{th=}")
        mu = mumu
    else:
        th = pl.strict_theta(~pl.Permutation(v))
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
            sumval = vpathsums.get(v, 0)
            if sumval == 0:
                continue
            for v2, vdiff, s in vpathdicts[index][v]:
                newpathsums[v2] = newpathsums.get(
                    v2,
                    0,
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
    return vpathsums.get(vmu, 0)

def schubpoly_classical_from_elems(v, var_x=None, var_y=None, elem_func=None):
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
            sumval = vpathsums.get(v, 0)
            if sumval == 0:
                continue
            for v2, vdiff, s in vpathdicts[index][v]:
                newpathsums[v2] = newpathsums.get(
                    v2,
                    0,
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
    return vpathsums.get(vmu, 0)


def schubpoly(v, var2=GeneratingSet("y"), var3=GeneratingSet("z"), start_var=1):
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


def skew_div_diff(u, w, poly):
    d = -1
    for i in range(len(w) - 1):
        if w[i] > w[i + 1]:
            d = i
            break
    d2 = -1
    for i in range(len(u) - 1):
        if u[i] > u[i + 1]:
            d2 = i
            break
    if d == -1:
        if d2 == -1:
            return poly
        return 0
    w2 = w.swap(d, d + 1)
    if d < len(u) - 1 and u[d] > u[d + 1]:
        u2 = u.swap(d, d + 1)
        return skew_div_diff(u2, w2, perm_act(poly, d + 1))
    return skew_div_diff(u, w2, div_diff(d + 1, poly))

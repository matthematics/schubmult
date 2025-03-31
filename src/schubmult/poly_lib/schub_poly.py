import sympy

from schubmult.schub_lib import pull_out_var


def perm_act(val, i, var2=None):
    subsdict = {var2[i]: var2[i + 1], var2[i + 1]: var2[i]}
    return sympy.ympify(val).subs(subsdict)


def div_diff(i, poly, var2=None):
    return sympy.sympify(
        sympy.div(sympy.sympify(poly - perm_act(poly, i)), sympy.sympify(var2[i] - var2[i + 1]))[0],
    )


def schubpoly(v, var2=None, var3=None, start_var=1):
    n = 0
    for j in range(len(v) - 2, -1, -1):
        if v[j] > v[j + 1]:
            n = j + 1
            break
    if n == 0:
        return 1
    lst = pull_out_var(n, v)
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

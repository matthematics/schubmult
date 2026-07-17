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


def _groth_plus(x1, y1, beta):
    return x1 + y1 + beta * x1 * y1


def _groth_div_diff(val, index, x, beta):
    from schubmult.rings.schubert.schubert_ring import DoubleSchubertElement, SingleSchubertRing

    ring = SingleSchubertRing(x)
    if not isinstance(val, DoubleSchubertElement):
        val = ring.from_expr(val)
    up_val = (S.One + beta * x[index + 1]) * val
    rval = ring.from_dict({w.swap(index - 1, index): coeff for w, coeff in up_val.items() if w[index - 1] > w[index]})
    return rval


@cache
def grothendieck_poly(perm, x, y, beta, keep_as_schub=False):
    from schubmult.combinatorics.permutation import Permutation

    if perm.inv == 0:
        return S.One
    n = len(perm)
    w0 = Permutation.w0(n)
    if perm == w0:
        return prod([_groth_plus(x[i], y[j], beta) for i in range(1, n) for j in range(1, n + 1 - i)])
    desc = min([i for i in range(n) if i not in (perm.descents())])
    result = _groth_div_diff(grothendieck_poly(perm.swap(desc, desc + 1), x, y, beta, keep_as_schub=True), desc + 1, x, beta)
    if keep_as_schub:
        return result
    return result.as_polynomial()


def to_groth(val, x, y, beta):
    from schubmult import uncode
    from schubmult.rings.polynomial_algebra import MonomialBasis, PolynomialAlgebra
    from schubmult.symbolic import expand

    if expand(val, deep=True) == S.Zero:
        return {}
    PA = PolynomialAlgebra(MonomialBasis(x))
    as_pol = PA.from_dict({k: v for k, v in PA.from_expr(val).items() if expand(v) != S.Zero})
    comp = min(as_pol.keys(), key=lambda a: (-sum(a), a))
    # ring = SingleSchubertRing(x)
    # schub_dict = ring.from_dict({k: v for k, v in ring.from_expr(val).items() if expand(v) != S.Zero})
    # min_perm = min(schub_dict.keys(), key=lambda a: (a.inv, a.pad_code(max(len(p) for p in schub_dict.keys()))))
    new_val = val - as_pol[comp] * grothendieck_poly(uncode(comp), x, y, beta, keep_as_schub=False)
    assert expand(PA.from_expr(new_val, length=len(comp)).get(comp, S.Zero), deep=True) == S.Zero, f"Expected zero but got {expand(new_val, deep=True)}"
    old = {k: v for k, v in to_groth(new_val, x, y, beta).items() if expand(v) != S.Zero}
    perm = uncode(comp)
    old[perm] = old.get(perm, S.Zero) + as_pol[comp]
    if old[perm] == S.Zero:
        del old[perm]
    return old


def groth_dict_to_poly(groth_dict, x, zz, beta):
    ret = S.Zero
    for perm, coeff in groth_dict.items():
        ret += coeff * grothendieck_poly(perm, x, zz, beta)
    return ret

@cache
def schub_elem_to_groth_elem_dict(the_perm, beta):
    from schubmult import RCGraph
    from schubmult.combinatorics.pipe_dream import PipeDream
    w0 = pl.Permutation.w0(len(the_perm))
    dct = {}
    for bpd in RCGraph.all_rc_graphs(the_perm, len(the_perm)):
        permo = PipeDream.from_rc_graph(bpd).co_pipe_dream().perm * w0
        dct[(permo.inv, permo.max_descent)] = dct.get((permo.inv, permo.max_descent), S.Zero) + (-beta) ** (permo.inv - the_perm.inv)
    return dct

@cache
def schub_elem_sym_to_groth_elem_sym_dict(p, k, beta):
    from schubmult import RCGraph, uncode
    from schubmult.combinatorics.pipe_dream import PipeDream

    the_perm = uncode([0] * (k - p) + [1] * p)
    dct = {}
    w0 = pl.Permutation.w0(len(the_perm))
    for bpd in RCGraph.all_rc_graphs(the_perm, len(the_perm)):
        permo = PipeDream.from_rc_graph(bpd).co_pipe_dream().perm * w0
        dct[(permo.inv, permo.max_descent)] = dct.get((permo.inv, permo.max_descent), S.Zero) + (-beta) ** (permo.inv - the_perm.inv)
    ret = dct
    return ret


# @cache
def _strip_isobaric(index, length, genset, beta, elem, backwards=False):
    from schubmult import uncode
    from schubmult.abc import E
    from schubmult.rings.schubert.nil_hecke import NilHeckeRing
    from schubmult.rings.schubert.schubert_ring import SingleSchubertRing

    ring = SingleSchubertRing(genset)
    nh = NilHeckeRing(genset)
    if backwards:
        operator = nh(~uncode([0] * (index - 1) + [length]))
    else:
        operator = nh(uncode([0] * (index - 1) + [length]))
    schub = ring.from_dict({k: v * beta ** (k.inv) for k, v in (elem * E(length, length, genset[index + 1 :], [S.NegativeOne]) * ring.one).items()})
    # schub = ring.from_expr(poly)
    return operator.apply(schub)

def isobar_it(i, genset, elem):
    from schubmult import Permutation
    from schubmult.rings.schubert.nil_hecke import NilHeckeRing
    from schubmult.rings.schubert.schubert_ring import SingleSchubertRing
    ring = SingleSchubertRing(genset)
    nh = NilHeckeRing(genset)

    operator = nh(Permutation.ref_product(i))
    schub = ring.from_expr((1 + genset[i + 1]) * genset[i]) * elem
    return operator.apply(schub)

def lascoux_poly(composition, genset):
    from .. import expand
    return expand(_lascoux_poly(tuple(composition), genset))

@cache
def _lascoux_poly(composition, genset):
    if all(composition[i] >= composition[i + 1] for i in range(len(composition) - 1)):
        return prod([genset[i + 1] ** composition[i] for i in range(len(composition)) if composition[i] > 0])
    for i in range(len(composition) - 1):
        if composition[i] < composition[i + 1]:
            return isobar_it(i + 1, genset, lascoux_poly(composition[:i] + (composition[i + 1], composition[i]) + composition[i + 2 :], genset))
    raise ValueError(f"Unexpected composition: {composition}")

@cache
def groth_elem_as_schub_dict(perm, beta):
    """Expand a Grothendieck basis element G_perm into Schubert basis terms.

    Uses the strip-isobaric construction from the groth_elem_as_schub method.
    Returns ``{Permutation: coeff}``.
    """
    from schubmult import WCGraph

    # perm = pl.Permutation(perm)
    # ring = SingleSchubertRing(genset)
    if perm.inv == 0:
        return {pl.Permutation([]): S.One}

    # w0 = pl.Permutation.w0(len(perm))
    # schub_elem = ring(w0)
    # bacon = ((~perm) * w0).trimcode
    # for i in range(len(bacon) - 1, -1, -1):
    #     if bacon[i] != 0:
    #         schub_elem = _strip_isobaric(i + 1, bacon[i], genset, beta, schub_elem, backwards=False)
    # return {k: v for k, v in schub_elem.items() if v != S.Zero}
    return WCGraph.groth_to_schub(perm, beta)


# (1+beta y)E(x,y/(1+beta y))
# yy
# beta/(1-yy)E(x,yy)


def groth_mul_full(perm_dict, p2, _x, _zz, beta):
    #from schubmult import Gx, Sx
    p2 = pl.Permutation(p2)
    # schub_elem2 = Sx.from_dict(groth_elem_as_schub_dict(p2, beta))
    # schub_elem1 = sum([v * Sx.from_dict(groth_elem_as_schub_dict(k, beta)) for k, v in perm_dict.items()])
    # schub_dict_to_groth_dict currently takes (base_groth, schub_dict, beta).
    # Keep groth_mul_full signature stable for existing call sites, but route
    # through the updated helper interface.
    return schub_dict_to_groth_dict(perm_dict, groth_elem_as_schub_dict(p2, beta), beta)
    # result = schub_elem1 * schub_elem2
    # return sum([v * Gx.from_dict(schub_elem_to_groth_elem_dict(k, beta)) for k, v in result.items()])

def schub_dict_to_groth_dict(base_groth, schub_dict, beta):
    # schub_elem_sym_as_groth_elem_sym_dict

    import sympy

    from schubmult import FactorialElemSym, Sx
    from schubmult.utils.perm_utils import add_perm_dict
    from schubmult.utils.schub_lib import groth_pieri_mul

    # schub_poly_basis = SchubertPolyBasis(x)
    # the_elem_dict = schub_poly_basis.transition_elementary({(perm, max()): coeff for perm, coeff in schub_dict.items()}, ElemSymPolyBasis(x))
    schub_elem_expr = Sx.from_dict(schub_dict).in_CEM_basis()

    def _mul_scalar(term_dict, scalar):
        if scalar == S.Zero:
            return {}
        return {perm: coeff * scalar for perm, coeff in term_dict.items()}

    def _apply_factorial_elem_sym(term_dict, degree, numvars):
        if degree == 0:
            return term_dict
        dctt = schub_elem_sym_to_groth_elem_sym_dict(degree, numvars, beta)
        build = {}
        for pair, coeff3 in dctt.items():
            pieri_piece = groth_pieri_mul(term_dict, *pair, beta)
            build = add_perm_dict(build, {perm: coeff3 * coeff for perm, coeff in pieri_piece.items()})
        return build

    def _flatten_factors(expr):
        """Return (scalar, [(deg, numvars), ...]) for multiplicative terms.

        Returns None when expression contains additive structure that would require
        distribution; caller may then use recursive fallback.
        """
        expr_sym = sympy.sympify(expr)
        if isinstance(expr_sym, FactorialElemSym):
            return S.One, [(expr_sym.degree, expr_sym.numvars)]

        if isinstance(expr, Add):
            return None

        if isinstance(expr, Mul):
            scalar = S.One
            factors = []
            for arg in expr.args:
                flattened = _flatten_factors(arg)
                if flattened is None:
                    return None
                arg_scalar, arg_factors = flattened
                scalar *= arg_scalar
                factors.extend(arg_factors)
            return scalar, factors

        if isinstance(expr, Pow):
            base, exponent = expr.args
            exponent_int = int(exponent)
            if exponent_int < 0 or exponent_int != exponent:
                raise ValueError(f"Unsupported exponent in CEM expression: {exponent}")
            if exponent_int == 0:
                return S.One, []
            flattened_base = _flatten_factors(base)
            if flattened_base is None:
                return None
            base_scalar, base_factors = flattened_base
            return base_scalar**exponent_int, base_factors * exponent_int

        return expr, []

    def _eval_expr(term_dict, expr):
        expr_sym = sympy.sympify(expr)

        if isinstance(expr_sym, FactorialElemSym):
            return _apply_factorial_elem_sym(term_dict, expr_sym.degree, expr_sym.numvars)

        if isinstance(expr, Add):
            out = {}
            for arg in expr.args:
                out = add_perm_dict(out, _eval_expr(term_dict, arg))
            return out

        if isinstance(expr, Mul):
            out = term_dict
            for arg in expr.args:
                out = _eval_expr(out, arg)
            return out

        if isinstance(expr, Pow):
            base, exponent = expr.args
            exponent_int = int(exponent)
            if exponent_int < 0 or exponent_int != exponent:
                raise ValueError(f"Unsupported exponent in CEM expression: {exponent}")
            out = term_dict
            for _ in range(exponent_int):
                out = _eval_expr(out, base)
            return out

        return _mul_scalar(term_dict, expr)

    # Fast path: compile top-level additive pieces into multiplicative factor lists.
    top_terms = schub_elem_expr.args if isinstance(schub_elem_expr, Add) else (schub_elem_expr,)
    compiled_terms = []
    fallback_terms = []
    for term_expr in top_terms:
        flattened = _flatten_factors(term_expr)
        if flattened is None:
            fallback_terms.append(term_expr)
        else:
            compiled_terms.append(flattened)

    ret = {}
    for p1, coeff0 in base_groth.items():
        subtotal = {}

        for scalar, factors in compiled_terms:
            term_dict = {p1: coeff0}
            if scalar != S.One:
                term_dict = _mul_scalar(term_dict, scalar)
            for degree, numvars in factors:
                term_dict = _apply_factorial_elem_sym(term_dict, degree, numvars)
            subtotal = add_perm_dict(subtotal, term_dict)

        for term_expr in fallback_terms:
            subtotal = add_perm_dict(subtotal, _eval_expr({p1: coeff0}, term_expr))

        ret = add_perm_dict(ret, subtotal)
    return ret


if __name__ == "__main__":
    from symengine import Symbol, expand

    from schubmult import Permutation, Sx, uncode
    from schubmult.abc import x
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet

    zz = ZeroGeneratingSet()
    # beta = Symbol("beta")
    beta = S.NegativeOne
    Permutation.print_as_code = True
    # print(to_groth(grothendieck_poly(Permutation([2, 3, 4, 1]), x, y, beta), x, y, beta))
    test_poly = Sx(uncode([0, 0, 1, 1, 1])).as_polynomial()
    grothy = to_groth(test_poly, x, zz, beta)
    fat_poly = sum([v * grothendieck_poly(k, x, zz, beta) for k, v in grothy.items()])
    assert expand(fat_poly - test_poly, deep=True) == S.Zero, f"Expected zero but got {expand(fat_poly - test_poly, deep=True)}"
    # print(to_groth(Sx(Permutation([2, 5, 1, 4, 3])).as_polynomial(), x, zz, S.One))

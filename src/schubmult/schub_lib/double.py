from bisect import bisect_left
from functools import cache, cached_property

import numpy as np
import psutil
import pulp as pu
import sympy
from cachetools import cached
from cachetools.keys import hashkey
from sortedcontainers import SortedList
from symengine import Add, Integer, Mul, Pow, sympify

from schubmult.perm_lib import (
    Permutation,
    code,
    cycle,
    dominates,
    inv,
    one_dominates,
    phi1,
    theta,
    uncode,
)
from schubmult.poly_lib.poly_lib import efficient_subs, elem_sym_func, elem_sym_poly, expand
from schubmult.poly_lib.schub_poly import schubpoly
from schubmult.poly_lib.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
from schubmult.schub_lib.schub_lib import (
    compute_vpathdicts,
    divdiffable,
    elem_sym_perms,
    elem_sym_perms_op,
    is_coeff_irreducible,
    is_split_two,
    pull_out_var,
    reduce_coeff,
    reduce_descents,
    try_reduce_u,
    try_reduce_v,
    will_formula_work,
)
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

zero = sympify(0)

logger = get_logger(__name__)


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



_vars = _gvars()


def count_sorted(mn, tp):
    index = bisect_left(mn, tp)
    ct = 0
    if mn[index] == tp:
        while index < len(mn) and mn[index] == tp:
            ct += 1
    return ct


# def E(p, k, varl=_vars.var2[1:]):
#     return elem_sym_poly(p, k, _vars.var1[1:], varl)


def single_variable(coeff_dict, varnum, var2=None):
    ret = {}
    for u in coeff_dict:
        if varnum - 1 < len(u):
            ret[u] = ret.get(u, 0) + var2[u[varnum - 1]] * coeff_dict[u]
        else:
            ret[u] = ret.get(u, 0) + var2[varnum] * coeff_dict[u]
        new_perms_k = elem_sym_perms(u, 1, varnum)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = elem_sym_perms(u, 1, varnum - 1)
        for perm, udiff in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u]
        for perm, udiff in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u]
    return ret


def single_variable_down(coeff_dict, varnum, var2=_vars.var2):
    ret = {}
    for u in coeff_dict:
        if varnum - 1 < len(u):
            ret[u] = ret.get(u, 0) + var2[u[varnum - 1]] * coeff_dict[u]
        else:
            ret[u] = ret.get(u, 0) + var2[varnum] * coeff_dict[u]
        new_perms_k = elem_sym_perms_op(u, 1, varnum)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = elem_sym_perms_op(u, 1, varnum - 1)
        for perm, udiff in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u]
        for perm, udiff in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u]
    return ret


def mult_poly_double(coeff_dict, poly, var_x=_vars.var1, var_y=_vars.var2):
    # try:
    #     poly = sympify(poly)
    # except SympifyError:
    #     poly = sympy.sympify(poly)
    #     var_x = tuple([sympy.sympify(v) for v in var_x])
    #     var_y = tuple([sympy.sympify(v) for v in var_y])
    #     return mult_poly_sympy(coeff_dict, poly, var_x=_vars.var1, var_y=_vars.var2)
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    if var_x.index(poly) != -1:
        return single_variable(coeff_dict, var_x.index(poly), var_y)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_double(ret, a, var_x, var_y)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_double(ret, base, var_x, var_y)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_double(coeff_dict, a, var_x, var_y))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


# def mult_poly_symy(coeff_dict, poly, var_x=_vars.sympy_var1, var_y=_vars.sympy_var2):


def mult_poly_down(coeff_dict, poly):
    if poly in _vars.var1:
        return single_variable_down(coeff_dict, _vars.var1.index(poly))
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_down(ret, a)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_down(ret, base)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_down(coeff_dict, a))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def nilhecke_mult(coeff_dict1, coeff_dict2):
    ret = {}
    for w in coeff_dict2:
        w1 = w
        inv_w1 = inv(w1)
        poly = coeff_dict2[w]
        did_mul = mult_poly_down(coeff_dict1, poly)
        for v in did_mul:
            v1 = [*v]
            addperm = v1 * w1
            if inv(addperm) == inv(v1) + inv_w1:
                toadd = addperm
                ret[toadd] = ret.get(toadd, 0) + did_mul[v]
    return ret


def forwardcoeff(u, v, perm, var2=None, var3=None):
    th = theta(v)
    muv = uncode(th)
    vmun1 = (~v) * muv

    w = perm * vmun1
    if inv(w) == inv(vmun1) + inv(perm):
        coeff_dict = schubmult_double_pair(u, muv, var2, var3)
        # logger.debug(f"{coeff_dict.get(w,0)=} {w=} {perm=} {vmun1=} {v=} {muv=}")
        return coeff_dict.get(w, 0)
    return 0


def dualcoeff(u, v, perm, var2=None, var3=None):
    if inv(u) == 0:
        # logger.debug("Recording line number")
        vp = v * (~perm)
        if inv(vp) == inv(v) - inv(perm):
            return schubpoly(vp, var2, var3)
    dpret = []
    ret = 0
    if dominates(u, perm):
        dpret = dualpieri(u, v, perm)
    else:
        # logger.debug("Recording line number")
        dpret = []
        # logger.debug("Recording line number")
        th = theta(u)
        muu = uncode(th)
        umun1 = (~u) * muu
        w = perm * umun1
        # logger.debug("spiggle")
        # logger.debug(f"{u=} {muu=} {v=} {w=} {perm=}")
        # logger.debug(f"{w=} {perm=}")
        if inv(w) == inv(umun1) + inv(perm):
            dpret = dualpieri(muu, v, w)
            # logger.debug(f"{muu=} {v=} {w=}")
            # logger.debug(f"{dpret=}")
    for vlist, vp in dpret:
        # logger.debug("Recording line number")
        toadd = 1
        for i in range(len(vlist)):
            for j in range(len(vlist[i])):
                toadd *= var2[i + 1] - var3[vlist[i][j]]
        toadd *= schubpoly(vp, var2, var3, len(vlist) + 1)
        ret += toadd
    return ret
    # logger.debug("Recording line number")
    # schub_val = schubmult_one(u, v, var2, var3)
    # val_ret = schub_val.get(perm, 0)
    # if expand(val - val_ret) != 0:
    #     # logger.debug(f"{schub_val=}")
    #     # logger.debug(f"{val=} {u=} {v=} {var2[1]=} {var3[1]=}  {perm=} {schub_val.get(perm,0)=}")
    # logger.debug(f"good to go {ret=}")


def dualpieri(mu, v, w):
    # logger.debug(f"dualpieri {mu=} {v=} {w=}")
    lm = code(~mu)
    cn1w = code(~w)
    while len(lm) > 0 and lm[-1] == 0:
        lm.pop()
    while len(cn1w) > 0 and cn1w[-1] == 0:
        cn1w.pop()
    if len(cn1w) < len(lm):
        return []
    for i in range(len(lm)):
        if lm[i] > cn1w[i]:
            return []
    c = Permutation([1, 2])
    for i in range(len(lm), len(cn1w)):
        c = cycle(i - len(lm) + 1, cn1w[i]) * c
    # c = permtrim(c)
    # logger.debug("Recording line number")
    res = [[[], v]]
    # logger.debug(f"{v=} {type(v)=}")
    for i in range(len(lm)):
        # logger.debug(f"{res=}")
        res2 = []
        for vlist, vplist in res:
            vp = vplist
            vpl = divdiffable(vp, cycle(lm[i] + 1, cn1w[i] - lm[i]))
            # logger.debug(f"{vpl=} {type(vpl)=}")
            if len(vpl) == 0:
                continue
            vl = pull_out_var(lm[i] + 1, vpl)
            # logger.debug(f"{vl=}")
            for pw, vpl2 in vl:
                res2 += [[[*vlist, pw], vpl2]]
        res = res2
    if len(lm) == len(cn1w):
        return res
    res2 = []
    for vlist, vplist in res:
        vp = vplist
        vpl = divdiffable(vp, c)
        if len(vpl) == 0:
            continue
        res2 += [[vlist, vpl]]
    # logger.debug(f"{res2=}")
    return res2


dimen = 0
monom_to_vec = {}


@cache
def schubmult_double_pair(perm1, perm2, var2=None, var3=None):
    return schubmult_double({perm1: 1}, perm2, var2, var3)


@cache
def schubmult_double_pair_generic(perm1, perm2):
    return schubmult_double({perm1: 1}, perm2, _vars.var_g1, _vars.var_g2)


def schubmult_double(perm_dict, v, var2=None, var3=None):
    if isinstance(var2, str):
        var2 = GeneratingSet(var2)
    if isinstance(var3, str):
        var3 = GeneratingSet(var3)
    perm_dict = {Permutation(k): v for k, v in perm_dict.items()}
    v = Permutation(v)
    vn1 = ~v
    th = theta(vn1)
    if len(th) == 0:
        return perm_dict
    if th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu, True)
    for u, val in perm_dict.items():
        inv_u = inv(u)
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                inv_up = inv(up)
                newperms = elem_sym_perms(
                    up,
                    min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                    th[index],
                )
                for up2, udiff in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, zero)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                zero,
                            ) + s * sumval * elem_sym_func(
                                th[index],
                                index + 1,
                                up,
                                up2,
                                v,
                                v2,
                                udiff,
                                vdiff,
                                var2,
                                var3,
                            )
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({Permutation(ep): vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schubmult_down(perm_dict, v, var2=None, var3=None):
    vn1 = ~v
    th = theta(vn1)
    if th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    ret_dict = {}

    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu, True)
    for u, val in perm_dict.items():
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                newperms = elem_sym_perms_op(up, mx_th, th[index])
                for up2, udiff in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, zero)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                zero,
                            ) + s * sumval * elem_sym_func(
                                th[index],
                                index + 1,
                                up2,
                                up,
                                v,
                                v2,
                                udiff,
                                vdiff,
                                var2,
                                var3,
                            )
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def poly_to_vec(poly, vec0=None, var3=_vars.var3):
    poly = expand(sympify(poly).xreplace({var3[1]: 0}))

    dc = poly.as_coefficients_dict()

    if vec0 is None:
        init_basevec(dc)

    vec = {}
    for mn in dc:
        cf = dc[mn]
        if cf == 0:
            continue
        cf = abs(int(cf))
        try:
            index = monom_to_vec[mn]
        except KeyError:
            return None
        if vec0 is not None and vec0[index] < cf:
            return None
        vec[index] = cf
    return vec


def shiftsub(pol, var2=_vars.var2):
    subs_dict = {var2[i]: var2[i + 1] for i in range(99)}
    return efficient_subs(sympify(pol), subs_dict)


def shiftsubz(pol, var3=_vars.var3):
    subs_dict = {var3[i]: var3[i + 1] for i in range(99)}
    return efficient_subs(sympify(pol), subs_dict)


def init_basevec(dc):
    global dimen, monom_to_vec, base_vec  # noqa: PLW0603
    monom_to_vec = {}
    index = 0
    for mn in dc:
        if dc[mn] == 0:
            continue
        monom_to_vec[mn] = index
        index += 1
    dimen = index
    base_vec = [0 for i in range(dimen)]


def split_flat_term(arg):
    arg = expand(arg)
    ys = []
    zs = []
    for arg2 in arg.args:
        if str(arg2).find("y") != -1:
            if isinstance(arg2, Mul):
                for i in range(int(arg2.args[0])):
                    ys += [arg2.args[1]]
            else:
                ys += [arg2]
        elif isinstance(arg2, Mul):
            for i in range(abs(int(arg2.args[0]))):
                zs += [-arg2.args[1]]
        else:
            zs += [arg2]
    return ys, zs


def is_flat_term(term):
    if isinstance(term, Integer) or isinstance(term, int):
        return True
    dc = expand(term).as_coefficients_dict()
    for t in dc:
        if str(t).count("y") + str(t).count("z") > 1 or str(t).find("**") != -1:
            return False
    return True


def flatten_factors(term):
    found_one = False
    if is_flat_term(term):
        return term, False
    if isinstance(term, Pow):
        if is_flat_term(term.args[0]) and len(term.args[0].args) > 2:
            ys, zs = split_flat_term(term.args[0])
            terms = [1]
            for i in range(len(ys)):
                terms2 = []
                for j in range(len(term.args[1])):
                    for t in terms:
                        terms2 += [t * (ys[i] + zs[i])]
                terms = terms2
            return Add(*terms)
        if is_flat_term(term.args[0]):
            return term, False
        return flatten_factors(term.args[0]) ** term.args[1], True
    if isinstance(term, Mul):
        terms = [1]
        for arg in term.args:
            terms2 = []
            if isinstance(arg, Add) and not is_flat_term(expand(arg)):
                found_one = True
                for term3 in terms:
                    for arg2 in arg.args:
                        flat, found = flatten_factors(arg2)
                        terms2 += [term3 * flat]
            elif isinstance(arg, Add) and is_flat_term(arg) and len(arg.args) > 2:
                found_one = True
                ys, zs = split_flat_term(arg)
                for term3 in terms:
                    for i in range(len(ys)):
                        terms2 += [term3 * (ys[i] + zs[i])]
            else:
                flat, found = flatten_factors(arg)
                if found:
                    found_one = True
                for term3 in terms:
                    terms2 += [term3 * flat]
            terms = terms2
        if len(terms) == 1:
            term = terms[0]
        else:
            term = Add(*terms)
        return term, found_one
    if isinstance(term, Add):
        res = 0
        for arg in term.args:
            flat, found = flatten_factors(arg)
            if found:
                found_one = True
            res += flat
        return res, found_one
    return None


def fres(v):
    for s in v.free_symbols:
        return s
    return None


def split_mul(arg0, var2=None, var3=None):
    monoms = SortedList()

    var2s = {fres(var2[i]): i for i in range(len(var2))}
    var3s = {fres(var3[i]): i for i in range(len(var3))}
    # print(f"{type(arg0)=} {arg0=}")
    if isinstance(arg0, Pow):
        arg = arg0
        arg2 = expand(arg.args[0])
        yval = arg2.args[0]
        zval = arg2.args[1]
        if str(yval).find("z") != -1:
            yval, zval = zval, yval
        if str(zval).find("-") != -1:
            zval = -zval
        if str(yval).find("-") != -1:
            yval = -yval
        tup = (var2s[fres(yval)], var3s[fres(zval)])
        for i in range(int(arg0.args[1])):
            monoms += [tup]
    else:
        for arg in arg0.args:
            if is_flat_term(arg):
                if isinstance(arg, Integer) or isinstance(arg, int):
                    continue
                arg = expand(arg)
                if arg == 0:
                    break
                yval = arg.args[0]
                zval = arg.args[1]
                if str(yval).find("z") != -1:
                    yval, zval = zval, yval
                if str(zval).find("-") != -1:
                    zval = -zval
                if str(yval).find("-") != -1:
                    yval = -yval
                monoms += [(var2s[fres(yval)], var3s[fres(zval)])]
            elif isinstance(arg, Pow):
                arg2 = arg.args[0]
                yval = arg2.args[0]
                zval = arg2.args[1]
                if str(yval).find("z") != -1:
                    yval, zval = zval, yval
                if str(zval).find("-") != -1:
                    zval = -zval
                if str(yval).find("-") != -1:
                    yval = -yval
                tup = (var2s[fres(yval)], var3s[fres(zval)])
                for i in range(int(arg.args[1])):
                    monoms += [tup]
    return monoms


def split_monoms(pos_part, var2, var3):
    arrs = SortedList()
    if isinstance(pos_part, Add):
        for arg0 in pos_part.args:
            monoms = split_mul(arg0, var2, var3)
            arrs += [monoms]
    elif isinstance(pos_part, Mul) or isinstance(pos_part, Pow):
        monoms = split_mul(pos_part, var2, var3)
        arrs += [monoms]
    else:
        return [pos_part]
    return arrs


def is_negative(term):
    sign = 1
    if isinstance(term, Integer) or isinstance(term, int):
        return term < 0
    if isinstance(term, Mul):
        for arg in term.args:
            if isinstance(arg, Integer):
                sign *= arg
            elif isinstance(arg, Add):
                if str(arg).find("-y") != -1:
                    sign *= -1
            elif isinstance(arg, Pow):
                mulsign = 1
                if str(arg.args[0]).find("-y") != -1:
                    mulsign = -1
                sign *= mulsign**term.index
    elif isinstance(term, Pow):
        mulsign = 1
        if str(term.args[0]).find("-y") != -1:
            mulsign = -1
        sign *= mulsign**term.index
    return sign < 0


def find_base_vectors(monom_list, var2, var3, depth):
    size = 0
    mn_fullcount = {}
    # pairs_checked = set()
    monom_list = {tuple(mn) for mn in monom_list}
    ct = 0
    while ct < depth and size != len(monom_list):
        size = len(monom_list)
        monom_list2 = set(monom_list)
        additional_set2 = set()
        for mn in monom_list:
            mncount = mn_fullcount.get(mn, {})
            if mncount == {}:
                for tp in mn:
                    mncount[tp] = mncount.get(tp, 0) + 1
                mn_fullcount[mn] = mncount
            for mn2 in monom_list:
                mn2count = mn_fullcount.get(mn2, {})
                if mn2count == {}:
                    for tp in mn2:
                        mn2count[tp] = mn2count.get(tp, 0) + 1
                    mn_fullcount[mn2] = mn2count
                num_diff = 0
                for tp in mncount:
                    pt = mn2count.get(tp, 0) - mncount[tp]
                    num_diff += abs(pt)
                    if num_diff > 1:
                        break
                if num_diff == 1:
                    diff_term1 = None
                    diff_term2 = None
                    for tp in mn2count:
                        if mn2count[tp] > mncount.get(tp, 0):
                            diff_term2 = tp
                            break
                    for tp2 in mncount:
                        if mncount[tp2] > mn2count.get(tp2, 0):
                            diff_term1 = tp2
                            break
                    # print(f"{mn,mn2}")
                    if diff_term1 is None or diff_term2 is None:
                        raise Exception(f"{mn=} {mn2=}")
                    if diff_term2[1] == diff_term1[1]:
                        continue
                    new_term1 = (diff_term1[0], diff_term2[1])
                    new_term2 = (diff_term2[0], diff_term1[1])
                    # mn3 = [*mn]
                    # mn4 = list(mn2)
                    index = bisect_left(mn, diff_term1)
                    mn3 = list(mn[:index]) + list(mn[index + 1 :])
                    index = bisect_left(mn3, new_term1)
                    mn3_t = tuple(mn3[:index] + [new_term1] + mn3[index:])
                    index2 = bisect_left(mn2, diff_term2)
                    mn4 = list(mn2[:index2]) + list(mn2[index2 + 1 :])
                    index2 = bisect_left(mn4, new_term2)
                    mn4_t = tuple(mn4[:index2] + [new_term2] + mn4[index2:])
                    if mn3_t not in monom_list2:
                        additional_set2.add(mn3_t)
                    monom_list2.add(mn3_t)
                    if mn4_t not in monom_list2:
                        additional_set2.add(mn4_t)
                    monom_list2.add(mn4_t)
        monom_list = monom_list2
        ct += 1
    ret = []
    for mn in monom_list:
        if len(mn) != len(set(mn)):
            continue
        res = 1
        for tp in mn:
            res *= var2[tp[0]] - var3[tp[1]]
        ret += [res]
    return ret, monom_list


def compute_positive_rep(val, var2=GeneratingSet("y"), var3=GeneratingSet("z"), msg=False, do_pos_neg=True):
    do_pos_neg = False
    notint = False
    try:
        val2 = int(expand(val))
        # val2 = expand(val)
    except Exception:
        notint = True
    if notint:
        frees = val.free_symbols
        # logger.debug(f"{frees=}")
        # logger.debug(f"{[type(s) for s in frees]=}")
        varsimp2 = [m for m in frees if var2.index(m) != -1]
        varsimp3 = [m for m in frees if var3.index(m) != -1]
        varsimp2.sort(key=lambda k: var2.index(k))
        varsimp3.sort(key=lambda k: var3.index(k))
        # logger.debug(f"{varsimp2=}")
        # logger.debug(f"{varsimp3=}")
        var22 = [sympy.sympify(v) for v in varsimp2]
        var33 = [sympy.sympify(v) for v in varsimp3]
        # var22 = [sympy.sympify(m) for m in varsimp2]
        # var33 = [sympy.sympify(m) for m in varsimp3]
        n1 = len(varsimp2)

        # for i in range(len(varsimp2)):
        #     varsimp2[i] = var2[var2list.index(varsimp2[i])]
        # for i in range(len(varsimp3)):
        #     varsimp3[i] = var3[var3list.index(varsimp3[i])]

        base_vectors = []
        base_monoms = []
        vec = poly_to_vec(val, None)

        if do_pos_neg:
            smp = val
            flat, found_one = flatten_factors(smp)
            while found_one:
                flat, found_one = flatten_factors(flat, varsimp2, varsimp3)
            pos_part = 0
            neg_part = 0
            if isinstance(flat, Add) and not is_flat_term(flat):
                for arg in flat.args:
                    if expand(arg) == 0:
                        continue
                    if not is_negative(arg):
                        pos_part += arg
                    else:
                        neg_part -= arg
                if neg_part == 0:
                    # print("no neg")
                    return pos_part
            depth = 1

            mons = split_monoms(pos_part, varsimp2, varsimp3)
            mons = {tuple(mn) for mn in mons}
            mons2 = split_monoms(neg_part, varsimp2, varsimp3)
            mons2 = {tuple(mn2) for mn2 in mons2}

            # mons2 = split_monoms(neg_part)
            # for mn in mons2:
            # if mn not in mons:
            # mons.add(mn)
            # print(mons)
            status = 0
            size = len(mons)
            while status != 1:
                base_monoms, mons = find_base_vectors(mons, mons2, varsimp2, varsimp3, depth)
                if len(mons) == size:
                    raise ValueError("Found counterexample")

                size = len(mons)
                base_vectors = []
                bad = False
                bad_vectors = []
                for i in range(len(base_monoms)):
                    vec0 = poly_to_vec(base_monoms[i], vec)
                    if vec0 is not None:
                        base_vectors += [vec0]
                    else:
                        bad_vectors += [i]
                for j in range(len(bad_vectors) - 1, -1, -1):
                    base_monoms.pop(bad_vectors[j])

                vrs = [pu.LpVariable(name=f"a{i}", lowBound=0, cat="Integer") for i in range(len(base_vectors))]
                lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
                lp_prob += 0
                eqs = [*base_vec]
                for j in range(len(base_vectors)):
                    for i in base_vectors[j]:
                        bvi = base_vectors[j][i]
                        if bvi == 1:
                            eqs[i] += vrs[j]
                        else:
                            eqs[i] += bvi * vrs[j]
                for i in range(dimen):
                    try:
                        lp_prob += eqs[i] == vec[i]
                    except TypeError:
                        bad = True
                        break
                if bad:
                    continue
                try:
                    solver = pu.PULP_CBC_CMD(msg=msg)
                    status = lp_prob.solve(solver)
                except KeyboardInterrupt:
                    current_process = psutil.Process()
                    children = current_process.children(recursive=True)
                    for child in children:
                        child_process = psutil.Process(child.pid)
                        child_process.terminate()
                        child_process.kill()
                    raise KeyboardInterrupt()
                status = lp_prob.status
        else:
            # logger.debug("this")
            val_poly = sympy.poly(expand(val), *var22, *var33)
            vec = poly_to_vec(val)
            mn = val_poly.monoms()
            L1 = tuple([0 for i in range(n1)])
            mn1L = []
            lookup = {}
            # logger.debug("this")
            for mm0 in mn:
                key = mm0[n1:]
                if key not in lookup:
                    lookup[key] = []
                mm0n1 = mm0[:n1]
                st = set(mm0n1)
                if len(st.intersection({0, 1})) == len(st) and 1 in st:
                    lookup[key] += [mm0]
                if mm0n1 == L1:
                    mn1L += [mm0]
            # logger.debug("this")
            for mn1 in mn1L:
                comblistmn1 = [1]
                for i in range(n1, len(mn1)):
                    if mn1[i] != 0:
                        arr = np.array(comblistmn1)
                        comblistmn12 = []
                        mn1_2 = (*mn1[n1:i], 0, *mn1[i + 1 :])
                        for mm0 in lookup[mn1_2]:
                            comblistmn12 += (
                                arr
                                * np.prod(
                                    [varsimp2[k] - varsimp3[i - n1] for k in range(n1) if mm0[k] == 1],
                                )
                            ).tolist()
                        comblistmn1 = comblistmn12
                for i in range(len(comblistmn1)):
                    b1 = comblistmn1[i]
                    vec0 = poly_to_vec(b1, vec)
                    if vec0 is not None:
                        base_vectors += [vec0]
                        base_monoms += [b1]
            vrs = [pu.LpVariable(name=f"a{i}", lowBound=0, cat="Integer") for i in range(len(base_vectors))]
            lp_prob = pu.LpProblem("Problem", pu.LpMinimize)
            lp_prob += 0
            eqs = [*base_vec]
            for j in range(len(base_vectors)):
                for i in base_vectors[j]:
                    bvi = base_vectors[j][i]
                    if bvi == 1:
                        eqs[i] += vrs[j]
                    else:
                        eqs[i] += bvi * vrs[j]
            for i in range(dimen):
                lp_prob += eqs[i] == vec[i]
            # logger.debug("I IS SOLVING")
            try:
                # logger.debug("I IS SOLVING BOLVING")
                solver = pu.PULP_CBC_CMD(msg=msg)
                status = lp_prob.solve(solver)
            except KeyboardInterrupt:
                current_process = psutil.Process()
                children = current_process.children(recursive=True)
                for child in children:
                    child_process = psutil.Process(child.pid)
                    child_process.terminate()
                    child_process.kill()
                raise KeyboardInterrupt()
        # print(f"{pos_part=}")
        # print(f"{neg_part=}")
        # else:
        # print(f"No dice {flat=}")
        # exit(1)
        # #val = pos_part - neg_part

        # depth+=1
        val2 = 0
        for k in range(len(base_vectors)):
            x = vrs[k].value()
            b1 = base_monoms[k]
            if x != 0 and x is not None:
                val2 += int(x) * b1
    # logger.debug(f"{val2=}")
    return val2


def posify_generic_partial(val, u2, v2, w2):
    val2 = val
    val = posify(val, u2, v2, w2, var2=_vars.var_g1, var3=_vars.var_g2, msg=True, do_pos_neg=False, sign_only=False, optimize=False)
    if expand(val - val2) != 0:
        # logger.debug("Warning, failed on a case")
        raise Exception(f"{val=} {val2=} {u2=} {v2=} {w2=}")
    # print("FROFL")
    return val


@cache
def schubmult_generic_partial_posify(u2, v2):
    return {w2: posify_generic_partial(val, u2, v2, w2) for w2, val in schubmult_double_pair_generic(u2, v2).items()}


@cached(
    cache={},
    key=lambda val, u2, v2, w2, var2=None, var3=None, msg=False, do_pos_neg=False, sign_only=False, optimize=True: hashkey(val, u2, v2, w2, var2, var3, msg, do_pos_neg, sign_only, optimize),
)
def posify(
    val,
    u2,
    v2,
    w2,
    var2=None,
    var3=None,
    msg=False,
    do_pos_neg=False,
    sign_only=False,
    optimize=True,
    n=_vars.n,
):
    # logger.debug(f"NEW {val=} {u2=} {v2=} {w2=}")
    oldval = val
    if inv(u2) + inv(v2) - inv(w2) == 0:
        # logger.debug(f"Hmm this is probably not or val inty true {val=}")
        return val

    if not sign_only and expand(val) == 0:
        # logger.debug(f"Hmm this is probably not true {u2=} {v2=} {w2=} {val=}")
        return 0
    # logger.debug("proceeding")
    u, v, w = u2, v2, w2
    # u, v, w = try_reduce_v(u2, v2, w2)
    if is_coeff_irreducible(u2, v2, w2):
        u, v, w = try_reduce_u(u2, v2, w2)
        if is_coeff_irreducible(u, v, w):
            u, v, w = u2, v2, w2
            if is_coeff_irreducible(u, v, w):
                w0 = w
                u, v, w = reduce_descents(u, v, w)
                if is_coeff_irreducible(u, v, w):
                    u, v, w = reduce_coeff(u, v, w)
                    if is_coeff_irreducible(u, v, w):
                        while is_coeff_irreducible(u, v, w) and w0 != w:
                            w0 = w
                            u, v, w = reduce_descents(u, v, w)
                            if is_coeff_irreducible(u, v, w):
                                u, v, w = reduce_coeff(u, v, w)

    if w != w2 and sign_only:
        # logger.debug(f"Return 0 ")
        return 0
    # logger.debug(f"Reduced to {u2=} {v2=} {w2=} {val=}")
    if is_coeff_irreducible(u, v, w):
        u3, v3, w3 = try_reduce_v(u, v, w)
        if not is_coeff_irreducible(u3, v3, w3):
            u, v, w = u3, v3, w3
        else:
            u3, v3, w3 = try_reduce_u(u, v, w)
            if not is_coeff_irreducible(u3, v3, w3):
                u, v, w = u3, v3, w3
    split_two_b, split_two = is_split_two(u, v, w)
    # logger.debug("Recording line number")
    if len([i for i in code(v) if i != 0]) == 1:
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        cv = code(v)
        for i in range(len(cv)):
            if cv[i] != 0:
                k = i + 1
                p = cv[i]
                break
        inv_u = inv(u)
        r = inv(w) - inv_u
        val = 0
        w2 = w
        hvarset = [w2[i] for i in range(min(len(w2), k))] + [i + 1 for i in range(len(w2), k)] + [w2[b] for b in range(k, len(u)) if u[b] != w2[b]] + [w2[b] for b in range(len(u), len(w2))]
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        return elem_sym_poly(
            p - r,
            k + p - 1,
            [-var3[i] for i in range(1, n)],
            [-var2[i] for i in hvarset],
        )
        # if expand(val - oldval) != 0:
        #     # logger.debug("This is bad")
        #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
    if will_formula_work(v, u) or dominates(u, w):
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        return dualcoeff(u, v, w, var2, var3)
        # if expand(val - oldval) != 0:
        # logger.debug("This is bad")
        # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=} {will_formula_work(v,u)=} {dominates(u,w)=}")
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
    if inv(w) - inv(u) == 1:
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        a, b = -1, -1
        for i in range(len(w)):
            if a == -1 and u[i] != w[i]:
                a = i
            elif (i >= len(u) and w[i] != i + 1) or (b == -1 and u[i] != w[i]):
                b = i
        arr = [[[], v]]
        d = -1
        for i in range(len(v) - 1):
            if v[i] > v[i + 1]:
                d = i + 1
        for i in range(d):
            arr2 = []
            if i in [a, b]:
                continue
            i2 = 1
            if i > b:
                i2 += 2
            elif i > a:
                i2 += 1
            for vr, v2 in arr:
                dpret = pull_out_var(i2, v2)
                for v3r, v3 in dpret:
                    arr2 += [[[*vr, v3r], v3]]
            arr = arr2
        val = 0
        for L in arr:
            v3 = L[-1]
            if v3[0] < v3[1]:
                continue
            v3 = v3.swap(0, 1)
            toadd = 1
            for i in range(d):
                if i in [a, b]:
                    continue
                i2 = i
                if i > b:
                    i2 = i - 2
                elif i > a:
                    i2 = i - 1
                oaf = L[0][i2]
                if i >= len(w):
                    yv = i + 1
                else:
                    yv = w[i]
                for j in range(len(oaf)):
                    toadd *= var2[yv] - var3[oaf[j]]
            toadd *= schubpoly(v3, [0, var2[w[a]], var2[w[b]]], var3)
            val += toadd
        # if expand(val - oldval) != 0:
        # logger.debug("This is bad")
        # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        # logger.debug(f"good to go {u2=} {v2=} {w2=}")
        return val
    # if split_two_b:
    #     # logger.debug("Recording line number")
    #     if sign_only:
    #         return 0
    #     cycles = split_two
    #     a1, b1 = cycles[0]
    #     a2, b2 = cycles[1]
    #     a1 -= 1
    #     b1 -= 1
    #     a2 -= 1
    #     b2 -= 1
    #     spo = sorted([a1, b1, a2, b2])
    #     real_a1 = min(spo.index(a1), spo.index(b1))
    #     real_a2 = min(spo.index(a2), spo.index(b2))
    #     real_b1 = max(spo.index(a1), spo.index(b1))
    #     real_b2 = max(spo.index(a2), spo.index(b2))

    #     good1 = False
    #     good2 = False
    #     if real_b1 - real_a1 == 1:
    #         good1 = True
    #     if real_b2 - real_a2 == 1:
    #         good2 = True
    #     a, b = -1, -1
    #     if good1 and not good2:
    #         a, b = min(a2, b2), max(a2, b2)
    #     if good2 and not good1:
    #         a, b = min(a1, b1), max(a1, b1)
    #     arr = [[[], v]]
    #     d = -1
    #     for i in range(len(v) - 1):
    #         if v[i] > v[i + 1]:
    #             d = i + 1
    #     for i in range(d):
    #         arr2 = []

    #         if i in [a1, b1, a2, b2]:
    #             continue
    #         i2 = 1
    #         i2 += len([aa for aa in [a1, b1, a2, b2] if i > aa])
    #         for vr, v2 in arr:
    #             dpret = pull_out_var(i2, v2)
    #             for v3r, v3 in dpret:
    #                 arr2 += [[[*vr, (v3r, i + 1)], v3]]
    #         arr = arr2
    #     val = 0

    #     if good1:
    #         arr2 = []
    #         for L in arr:
    #             v3 = L[-1]
    #             if v3[real_a1] < v3[real_b1]:
    #                 continue
    #             v3 = v3.swap(real_a1, real_b1)
    #             arr2 += [[L[0], v3]]
    #         arr = arr2
    #         if not good2:
    #             for i in range(4):
    #                 arr2 = []

    #                 if i in [real_a2, real_b2]:
    #                     continue
    #                 if i == real_a1:
    #                     var_index = min(a1, b1) + 1
    #                 elif i == real_b1:
    #                     var_index = max(a1, b1) + 1
    #                 i2 = 1
    #                 i2 += len([aa for aa in [real_a2, real_b2] if i > aa])
    #                 for vr, v2 in arr:
    #                     dpret = pull_out_var(i2, v2)
    #                     for v3r, v3 in dpret:
    #                         arr2 += [[[*vr, (v3r, var_index)], v3]]
    #                 arr = arr2
    #     if good2:
    #         arr2 = []
    #         for L in arr:
    #             v3 = L[-1]
    #             try:
    #                 if v3[real_a2] < v3[real_b2]:
    #                     continue
    #                 v3 = v3.swap(real_a2, real_b2)
    #             except IndexError:
    #                 continue
    #             arr2 += [[L[0], v3]]
    #         arr = arr2
    #         if not good1:
    #             for i in range(4):
    #                 arr2 = []

    #                 if i in [real_a1, real_b1]:
    #                     continue
    #                 i2 = 1
    #                 i2 += len([aa for aa in [real_a1, real_b1] if i > aa])
    #                 if i == real_a2:
    #                     var_index = min(a2, b2) + 1
    #                 elif i == real_b2:
    #                     var_index = max(a2, b2) + 1
    #                 for vr, v2 in arr:
    #                     dpret = pull_out_var(i2, v2)
    #                     for v3r, v3 in dpret:
    #                         arr2 += [[[*vr, (v3r, var_index)], v3]]
    #                 arr = arr2

    #         for L in arr:
    #             v3 = L[-1]
    #             tomul = 1
    #             doschubpoly = True
    #             if (not good1 or not good2) and v3[0] < v3[1] and (good1 or good2):
    #                 continue
    #             if (good1 or good2) and (not good1 or not good2):
    #                 v3 = v3.swap(0, 1)
    #             elif not good1 and not good2:
    #                 doschubpoly = False
    #                 if v3[0] < v3[1]:
    #                     dual_u = uncode([2, 0])
    #                     dual_w = Permutation([4, 2, 1, 3])
    #                     coeff = perm_act(dualcoeff(dual_u, v3, dual_w, var2, var3), 2, var2)

    #                 elif len(v3) < 3 or v3[1] < v3[2]:
    #                     if len(v3) <= 3 or v3[2] < v3[3]:
    #                         coeff = 0
    #                         continue
    #                     v3 = v3.swap(0, 1).swap(2, 3)
    #                     coeff = perm_act(schubpoly(v3, var2, var3), 2, var2)
    #                 elif len(v3) <= 3 or v3[2] < v3[3]:
    #                     if len(v3) <= 3:
    #                         v3 += [4]
    #                     v3 = v3.swap(2, 3)
    #                     coeff = perm_act(
    #                         posify(
    #                             schubmult_one(Permutation([1, 3, 2]), v3, var2, var3).get(
    #                                 Permutation([2, 4, 3, 1]),
    #                                 0,
    #                             ),
    #                             Permutation([1, 3, 2]),
    #                             v3,
    #                             Permutation([2, 4, 3, 1]),
    #                             var2,
    #                             var3,
    #                             msg,
    #                             do_pos_neg,
    #                             optimize=optimize,
    #                         ),
    #                         2,
    #                         var2,
    #                     )
    #                     # logger.debug(f"{coeff=}")
    #                 else:
    #                     coeff = perm_act(
    #                         schubmult_one(Permutation([1, 3, 2]), v3, var2, var3).get(
    #                             Permutation([2, 4, 1, 3]),
    #                             0,
    #                         ),
    #                         2,
    #                         var2,
    #                     )
    #                 # logger.debug(f"{coeff=}")
    #                 # if expand(coeff) == 0:
    #                 #     # logger.debug("coeff 0 oh no")
    #                 tomul = sympify(coeff)
    #             toadd = 1
    #             for i in range(len(L[0])):
    #                 var_index = L[0][i][1]
    #                 oaf = L[0][i][0]
    #                 if var_index - 1 >= len(w):
    #                     yv = var_index
    #                 else:
    #                     yv = w[var_index - 1]
    #                 for j in range(len(oaf)):
    #                     toadd *= var2[yv] - var3[oaf[j]]
    #             if (not good1 or not good2) and (good1 or good2):
    #                 varo = [0, var2[w[a]], var2[w[b]]]
    #             else:
    #                 varo = [0, *[var2[w[spo[k]]] for k in range(4)]]
    #             if doschubpoly:
    #                 toadd *= schubpoly(v3, varo, var3)
    #             else:
    #                 subs_dict3 = {var2[i]: varo[i] for i in range(len(varo))}
    #                 toadd *= efficient_subs(tomul, subs_dict3)
    #             val += toadd
    #             # logger.debug(f"accum {val=}")
    #         #logger.debug(f"{expand(val-oldval)=}")
    #         # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
    #         return val
    if will_formula_work(u, v):
        # logger.debug("Recording line number")
        if sign_only:
            return 0
        # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        return forwardcoeff(u, v, w, var2, var3)
        # if expand(val - oldval) != 0:
        #     # logger.debug("This is bad")
        #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
    # logger.debug("Recording line number")
    # c01 = code(u)
    # c02 = code(w)
    # c03 = code(v)

    c1 = code(~u)
    c2 = code(~w)

    if one_dominates(u, w):
        if sign_only:
            return 0
        while c1[0] != c2[0]:
            w = w.swap(c2[0] - 1, c2[0])
            v = v.swap(c2[0] - 1, c2[0])
            # w[c2[0] - 1], w[c2[0]] = w[c2[0]], w[c2[0] - 1]
            # v[c2[0] - 1], v[c2[0]] = v[c2[0]], v[c2[0] - 1]
            # w = tuple(w)
            # v = tuple(v)
            c2 = code(~w)
            # c03 = code(v)
            # c01 = code(u)
            # c02 = code(w)
        # if is_reducible(v):
        #     # logger.debug("Recording line number")
        #     if sign_only:
        #         return 0
        #     newc = []
        #     elemc = []
        #     for i in range(len(c03)):
        #         if c03[i] > 0:
        #             newc += [c03[i] - 1]
        #             elemc += [1]
        #         else:
        #             break
        #     v3 = uncode(newc)
        #     coeff_dict = schubmult_one(
        #         u,
        #         uncode(elemc),
        #         var2,
        #         var3,
        #     )
        #     val = 0
        #     for new_w in coeff_dict:
        #         tomul = coeff_dict[new_w]
        #         newval = schubmult_one(new_w, uncode(newc), var2, var3).get(
        #             w,
        #             0,
        #         )
        #         # logger.debug(f"Calling posify on {newval=} {new_w=} {uncode(newc)=} {w=}")
        #         newval = posify(newval, new_w, uncode(newc), w, var2, var3, msg, do_pos_neg, optimize=optimize)
        #         val += tomul * shiftsubz(newval)
        #     # if expand(val - oldval) != 0:
        #     #     # logger.debug("This is bad")
        #     #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        #         # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        #     return val
        # removed, iffy (hard to implement)
        # if c01[0] == c02[0] and c01[0] != 0:
        #     # logger.debug("Recording line number")
        #     if sign_only:
        #         return 0
        #     varl = c01[0]
        #     u3 = uncode([0] + c01[1:])
        #     w3 = uncode([0] + c02[1:])
        #     val = 0
        #     val = schubmult_one(u3, v, var2, var3).get(
        #         w3,
        #         0,
        #     )
        #     # logger.debug(f"Calling posify on {val=} {u3=} {v=} {w3=}")
        #     val = posify(val, u3, v, w3, var2, var3, msg, do_pos_neg, optimize=optimize)
        #     for i in range(varl):
        #         val = perm_act(val, i + 1, var2)
        #     # if expand(val - oldval) != 0:
        #     #     # logger.debug("This is bad")
        #     #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=}")
        #     # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
        #     return val
        if c1[0] == c2[0]:
            # logger.debug("Recording line number")
            if sign_only:
                return 0
            vp = pull_out_var(c1[0] + 1, v)
            u3 = phi1(u)
            w3 = phi1(w)
            val = 0
            for arr, v3 in vp:
                tomul = 1
                for i in range(len(arr)):
                    tomul *= var2[1] - var3[arr[i]]

                val2 = schubmult_double_pair(u3, v3, var2, var3).get(
                    w3,
                    0,
                )
                val2 = posify(val2, u3, v3, w3, var2, var3, msg, do_pos_neg, optimize=optimize)
                val += tomul * shiftsub(val2)
            # if expand(val - oldval) != 0:
            #     # logger.debug("This is bad")
            #     # logger.debug(f"{u2=} {v2=} {w2=} {val=} {oldval=")
            # logger.debug(f"Returning {u2=} {v2=} {w2=} {val=}")
            return val
    # logger.debug("Fell all the way through. Cleverness did not save us")
    if not sign_only:
        # logger.debug("Recording line number")
        if optimize:
            if inv(u) + inv(v) - inv(w) == 1:
                val2 = compute_positive_rep(val, var2, var3, msg, False)
            else:
                val2 = compute_positive_rep(val, var2, var3, msg, do_pos_neg)
            if val2 is not None:
                val = val2
            return val
        # logger.debug("RETURNINGOLDVAL")
        return oldval
    # logger.debug("Recording line number")
    d = expand(val).as_coefficients_dict()
    for v in d.values():
        if v < 0:
            return -1
    return 1

    #     if expand(val - oldval) !=0:
    #         # logger.debug("NONONOONOO")
    #         raise Exception
    # return val


def schub_coprod_double(mperm, indices, var2=_vars.var2, var3=_vars.var3):
    indices = sorted(indices)
    subs_dict_coprod = {}
    k = len(indices)
    n = len(mperm)
    kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    max_required = max([kcd[i] + i for i in range(len(kcd))])
    kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
    N = len(kcd)
    kperm = ~uncode(kcd2)
    inv_kperm = inv(kperm)
    vn = GeneratingSet("soible")

    for i in range(1, N * 2 + 1):
        if i <= N:
            subs_dict_coprod[vn[i]] = var2[i]
        else:
            subs_dict_coprod[vn[i]] = var3[i - N]

    coeff_dict = {kperm: 1}
    coeff_dict = schubmult_double(coeff_dict, mperm, vn, var2)

    inverse_kperm = ~kperm

    ret_dict = {}
    for perm in coeff_dict:
        downperm = perm * inverse_kperm
        if inv(downperm) == inv(perm) - inv_kperm:
            flag = True
            for i in range(N):
                if downperm[i] > N:
                    flag = False
                    break
            if not flag:
                continue
            firstperm = Permutation(downperm[0:N])
            secondperm = Permutation([downperm[i] - N for i in range(N, len(downperm))])

            val = efficient_subs(sympify(coeff_dict[perm]), subs_dict_coprod)

            key = (firstperm, secondperm)
            ret_dict[key] = val

    return ret_dict

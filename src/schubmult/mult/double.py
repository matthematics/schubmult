from bisect import bisect_left
from functools import cache

from schubmult.rings.poly_lib import _vars, efficient_subs, elem_sym_func, elem_sym_poly
from schubmult.rings.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
from schubmult.schub_lib.perm_lib import (
    Permutation,
    inv,
    theta,
    uncode,
)
from schubmult.schub_lib.schub_poly import elem_func_func_mul
from schubmult.symbolic import Add, Mul, Pow, S, expand, expand_func, sympify
from schubmult.symmetric_polynomials import FactorialElemSym
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict
from schubmult.utils.schub_lib import (
    compute_vpathdicts,
    elem_sym_perms,
    elem_sym_perms_op,
    elem_sym_positional_perms,
    pull_out_var,
)

zero = sympify(0)

logger = get_logger(__name__)


def count_sorted(mn, tp):
    index = bisect_left(mn, tp)
    ct = 0
    if mn[index] == tp:
        while index < len(mn) and mn[index] == tp:
            ct += 1
    return ct


# def E(p, k, varl=None):
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


def single_variable_down(coeff_dict, varnum, var2=None):
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


def mult_poly_double(coeff_dict, poly, var_x=None, var_y=None):
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


def mult_poly_double_alt(coeff_dict, poly, var_x=None, var_y=None):
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    if var_x.index(poly) != -1:
        return single_variable(coeff_dict, var_x.index(poly), var_y)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            s_d = mult_poly_double_alt({Permutation([]): S.One}, a, var_x, var_y)
            ret = schubmult_double_dict(ret, s_d, var_y, var_y)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        s_d = mult_poly_double_alt({Permutation([]): S.One}, base, var_x, var_y)
        for i in range(int(exponent)):
            ret = schubmult_double_dict(ret, s_d, var_y, var_y)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_double_alt(coeff_dict, a, var_x, var_y))
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


@cache
def schubmult_double_pair(perm1, perm2, var2=None, var3=None):
    return schubmult_double({perm1: S.One}, perm2, var2, var3)


@cache
def schubmult_double_pair_generic(perm1, perm2):
    return schubmult_double({perm1: S.One}, perm2, _vars.var_g1, _vars.var_g2)


@cache
def schubmult_double_pair_generic_alt(perm1, perm2):
    return {k: expand_func(expand(v)) for k,v in schubmult_double_alt_from_elems({perm1: S.One}, perm2, _vars.var_g1, _vars.var_g2, elem_func=FactorialElemSym).items()}


def schubmult_double_dict(perm_dict1, perm_dict2, var2=None, var3=None):
    ret = {}
    for k, v in perm_dict2.items():
        ret = add_perm_dict(ret, {k2: v2 * v for k2, v2 in schubmult_double(perm_dict1, k, var2, var3).items()})
    return ret


def schubmult_double(perm_dict, v, var2=None, var3=None):
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
    vpathdicts = compute_vpathdicts(th, vmu)
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


def schubmult_double_alt(perm_dict, v, var2=None, var3=None, index=1):
    if v.inv == 0:
        return perm_dict
        # ret = S.Zero
    ret_dict = {}
    L = pull_out_var(1, ~v)
    for index_list, new_v in L:
        interim_dict = {}
        for u, val in perm_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_sym_poly(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        ret_dict = add_perm_dict(ret_dict, schubmult_double_alt(interim_dict, ~new_v, var2, var3, index + 1))
    return ret_dict

# forwards backwards
def schubmult_double_alt_from_elems_forwards(perm_dict, v, var2=None, var3=None, index=1, elem_func=None):
    if v.inv == 0:
        return perm_dict
        # ret = S.Zero
    ret_dict = {}
    L = pull_out_var(1, ~v)
    for index_list, new_v in L:
        interim_dict = {}
        for u, val in perm_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        ret_dict = add_perm_dict(ret_dict, schubmult_double_alt_from_elems_forwards(interim_dict, ~new_v, var2, var3, index + 1, elem_func))
    return ret_dict

# backwards mul after
# def schubmult_double_alt_from_elems(perm_dict, v, var2=None, var3=None, elem_func=None):
#     if v.inv == 0:
#         return perm_dict
#     ret_dict = {}
#     index = max((~v).descents()) + 1
#     L = pull_out_var(index, ~v)
#     for index_list, new_v in L:
#         interim_dict = {}
#         for u, val in perm_dict.items():
#             new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
#             for new_perm, p, sgn in new_perms:
#                 interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
#                     len(index_list) - p,
#                     len(index_list) - p,
#                     [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
#                     [var3[index]],
#                 )
#         ret_dict = add_perm_dict(ret_dict, schubmult_double_alt_from_elems(interim_dict, ~new_v, var2, var3, elem_func))
#     return ret_dict

#backwards mul before
def schubmult_double_alt_from_elems_backwards(perm_dict, v, var2=None, var3=None, elem_func=None):
    if v.inv == 0:
        return perm_dict
    ret_dict = {}
    index = max((~v).descents()) + 1
    L = pull_out_var(index, ~v)
    _cache = {}
    for index_list, new_v in L:
        if new_v not in _cache:
            _cache[new_v] = schubmult_double_alt_from_elems_backwards(perm_dict, ~new_v, var2, var3, elem_func)
        start_dict = _cache[new_v]
        # start_dict = perm_dict
        interim_dict = {}
        for u, val in start_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        ret_dict = add_perm_dict(ret_dict, interim_dict)
        # ret_dict = add_perm_dict(ret_dict,schubmult_double_alt_from_elems_backwards(interim_dict, ~new_v, var2, var3, elem_func))
    return ret_dict

def schubmult_double_alt_from_elems_backwards_backwards(perm_dict, v, var2=None, var3=None, elem_func=None):
    if v.inv == 0:
        return perm_dict
    ret_dict = {}
    index = max((~v).descents()) + 1
    L = pull_out_var(index, ~v)
    #_cache = {}
    for index_list, new_v in L:
        # if new_v not in _cache:
        #     _cache[new_v] = schubmult_double_alt_from_elems_backwards(perm_dict, ~new_v, var2, var3, elem_func)
        start_dict = perm_dict
        # start_dict = perm_dict
        interim_dict = {}
        for u, val in start_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        # ret_dict = add_perm_dict(ret_dict, interim_dict)
        ret_dict = add_perm_dict(ret_dict,schubmult_double_alt_from_elems_backwards_backwards(interim_dict, ~new_v, var2, var3, elem_func))
    return ret_dict


schubmult_double_alt_from_elems = schubmult_double_alt_from_elems_backwards


def schubmult_double_from_elems(perm_dict, v, var2=None, var3=None, elem_func=None):
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
    vpathdicts = compute_vpathdicts(th, vmu)
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
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                0,
                            ) + s * sumval * elem_func_func_mul(
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
                                elem_func=elem_func,
                            )
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({Permutation(ep): vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schubmult_double_down(perm_dict, v, var2=None, var3=None):
    vn1 = ~v
    th = theta(vn1)
    if len(th) == 0 or th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    ret_dict = {}

    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
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


def schub_coprod_double(mperm, indices, var2=None, var3=None):
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

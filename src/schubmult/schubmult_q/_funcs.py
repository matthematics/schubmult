from functools import cached_property

from symengine import Add, Mul, Pow

from schubmult.perm_lib.perm_lib import (
    Permutation,
    add_perm_dict,
    code,
    count_bruhat,
    count_less_than,
    cycle,
    dominates,
    ensure_perms,
    get_cycles,
    getpermval,
    has_bruhat_ascent,
    has_bruhat_descent,
    inv,
    inverse,
    is_parabolic,
    longest_element,
    medium_theta,
    mu_A,
    mulperm,
    old_code,
    omega,
    one_dominates,
    p_trans,
    permtrim,
    permtrim_list,
    phi1,
    sg,
    split_perms,
    strict_theta,
    theta,
    trimcode,
    uncode,
)
from schubmult.poly_lib import GeneratingSet, is_indexed
from schubmult.schub_lib import (
    check_blocks,
    compute_vpathdicts,
    divdiffable,
    double_elem_sym_q,
    elem_sym_perms,
    elem_sym_perms_op,
    elem_sym_perms_q,
    elem_sym_perms_q_op,
    is_coeff_irreducible,
    is_hook,
    is_reducible,
    is_split_two,
    kdown_perms,
    pull_out_var,
    reduce_coeff,
    reduce_descents,
    reduce_q_coeff,
    reduce_q_coeff_u_only,
    try_reduce_u,
    try_reduce_v,
    will_formula_work,
)


class _gvars:
    @cached_property
    def n(self):
        return 100

    @cached_property
    def var_x(self):
        return GeneratingSet("x")

    @cached_property
    def q_var(self):
        return GeneratingSet("q")


_vars = _gvars()


def single_variable(coeff_dict, varnum, var_q=_vars.q_var):
    ret = {}
    for u in coeff_dict:
        new_perms_k = elem_sym_perms_q(u, 1, varnum, var_q)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = elem_sym_perms_q(u, 1, varnum - 1, var_q)
        for perm, udiff, mul_val in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u] * mul_val
        for perm, udiff, mul_val in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u] * mul_val
    return ret


def mult_poly(coeff_dict, poly, var_x=_vars.var_x, var_q=_vars.q_var):
    if is_indexed(poly) and poly.base==var_x:
        return single_variable(coeff_dict, poly.args[1], var_q=var_q)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly(ret, a, var_x, var_q=var_q)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly(ret, base, var_x, var_q=var_q)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly(coeff_dict, a, var_x, var_q=var_q))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def schubmult_db(perm_dict, v, q_var=_vars.q_var):
    if inv(v) == 0:
        return perm_dict
    th = medium_theta(~v)
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    mu = permtrim(uncode(th))
    vmu = permtrim(v*mu)
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}

    thL = len(th)
    # if thL!=2 and len(set(thL))!=1:
    # raise ValueError("Not what I can do")
    vpathdicts = compute_vpathdicts(th, vmu, True)
    # print(f"{vpathdicts=}")
    for u, val in perm_dict.items():
        inv_u = inv(u)
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            if index > 0 and th[index - 1] == th[index]:
                continue
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            if index < len(th) - 1 and th[index] == th[index + 1]:
                mx_th1 = 0
                for vp in vpathdicts[index + 1]:
                    for v2, vdiff, s in vpathdicts[index + 1][vp]:
                        mx_th1 = max(mx_th1, th[index + 1] - vdiff)
                newpathsums = {}
                for up in vpathsums:
                    newpathsums0 = {}
                    inv_up = inv(up)
                    newperms = double_elem_sym_q(up, mx_th, mx_th1, th[index], q_var)
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff2, s2 in vpathdicts[index][v]:
                            for up1, udiff1, mul_val1 in newperms:
                                if (up1, udiff1, mul_val1) not in newpathsums0:
                                    newpathsums0[(up1, udiff1, mul_val1)] = {}
                                if udiff1 + vdiff2 == th[index]:
                                    newpathsums0[(up1, udiff1, mul_val1)][v2] = (
                                        newpathsums0[(up1, udiff1, mul_val1)].get(
                                            v2,
                                            0,
                                        )
                                        + s2 * sumval * mul_val1
                                    )

                    for up1, udiff1, mul_val1 in newpathsums0:
                        for v in vpathdicts[index + 1]:
                            sumval = newpathsums0[(up1, udiff1, mul_val1)].get(v, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff2, s2 in vpathdicts[index + 1][v]:
                                for up2, udiff2, mul_val2 in newperms[(up1, udiff1, mul_val1)]:
                                    if up2 not in newpathsums:
                                        newpathsums[up2] = {}
                                    if udiff2 + vdiff2 == th[index + 1]:
                                        newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s2 * sumval * mul_val2
            else:
                newpathsums = {}
                for up in vpathsums:
                    inv_up = inv(up)
                    newperms = elem_sym_perms_q(
                        up,
                        min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                        th[index],
                        q_var,
                    )
                    for up2, udiff, mul_val in newperms:
                        if up2 not in newpathsums:
                            newpathsums[up2] = {}
                        for v in vpathdicts[index]:
                            sumval = vpathsums[up].get(v, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff, s in vpathdicts[index][v]:
                                if udiff + vdiff == th[index]:
                                    newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s * sumval * mul_val
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict(
            {ep: vpathsums[ep].get(toget, 0) for ep in vpathsums},
            ret_dict,
        )
    return ret_dict


def schubmult(perm_dict, v):
    th = strict_theta(inverse(v))
    mu = permtrim(uncode(th))
    vmu = permtrim(mulperm([*v], mu))
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu, True)
    for u, val in perm_dict.items():
        inv_u = inv(u)
        vpathsums = {u: {(1, 2): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                inv_up = inv(up)
                newperms = elem_sym_perms_q(
                    up,
                    min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                    th[index],
                )
                for up2, udiff, mul_val in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            if udiff + vdiff == th[index]:
                                newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s * sumval * mul_val
            vpathsums = newpathsums
        toget = tuple(vmu)
        ret_dict = add_perm_dict(
            {ep: vpathsums[ep].get(toget, 0) for ep in vpathsums},
            ret_dict,
        )
    return ret_dict


def grass_q_replace(perm, k, d, n):
    if k - d < 0:
        return None
    cd = code(perm)
    for i in range(k - d, k):
        if i >= len(cd) or cd[i] < d:
            return None
    grass_rep = [0 for i in range(n)]
    perm2 = [*perm] + [i + 1 for i in range(len(perm), n)]
    for i in range(k, n):
        grass_rep[perm2[i] - 1] = 2
    num_0 = 0
    # print(f"{grass_rep=} {d=}")
    for i in range(len(grass_rep) - 1, -1, -1):
        if num_0 == d:
            break
        if grass_rep[i] == 0:
            grass_rep[i] = 1
            num_0 += 1
    num_2 = 0
    for i in range(len(grass_rep)):
        if num_2 == d:
            break
        if grass_rep[i] == 2:
            grass_rep[i] = 1
            num_2 += 1
    # print(f"New {grass_rep=}")
    k1 = k - d
    k2 = k + d
    pos_1 = 0
    pos_2 = 0
    pos_3 = 0
    new_perm = [0 for i in range(n)]
    for i in range(len(grass_rep)):
        if grass_rep[i] == 0:
            new_perm[pos_1] = i + 1
            pos_1 += 1
        if grass_rep[i] == 1:
            new_perm[k1 + pos_2] = i + 1
            pos_2 += 1
        if grass_rep[i] == 2:
            new_perm[k2 + pos_3] = i + 1
            pos_3 += 1
    return tuple(permtrim(new_perm))


def to_two_step(perm, k1, k2, n):
    rep = [0 for i in range(n)]
    perm2 = [*perm] + [i + 1 for i in range(len(perm), n)]
    for i in range(n):
        if i < k1:
            rep[perm2[i] - 1] = 0
        elif i < k2:
            rep[perm2[i] - 1] = 1
        else:
            rep[perm2[i] - 1] = 2
    return rep

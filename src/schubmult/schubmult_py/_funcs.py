from functools import cached_property

# from symengine import Add, Mul, Pow, symarray
from symengine import Add, Mul, Pow

from schubmult.perm_lib import (
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
from schubmult.poly_lib import (
    GeneratingSet,
    call_zvars,
    elem_sym_func,
    elem_sym_func_q,
    elem_sym_poly,
    elem_sym_poly_q,
    is_indexed,
    q_vector,
)
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


_vars = _gvars()


def single_variable(coeff_dict, varnum):
    ret = {}
    for u in coeff_dict:
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


def mult_poly(coeff_dict, poly, var_x=_vars.var_x):
    if is_indexed(poly) and poly.base==var_x:
        return single_variable(coeff_dict, poly.args[0])
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly(ret, a, var_x)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly(ret, base, var_x)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly(coeff_dict, a, var_x))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def schubmult(perm_dict, v):
    v = Permutation(v)
    vn1 = ~v
    th = theta(vn1)
    if th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v*mu
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    vpathdicts = compute_vpathdicts(th, vmu)

    mx_th = [0 for i in range(len(th))]
    for index in range(len(th)):
        for vp in vpathdicts[index]:
            for v2, vdiff, s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    for u, val in perm_dict.items():
        inv_u = inv(u)
        vpathsums = {u: {Permutation([1,2]): val}}

        for index in range(len(th)):
            newpathsums = {}
            for up in vpathsums:
                inv_up = inv(up)
                newperms = elem_sym_perms(
                    up,
                    min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                    th[index],
                )
                for vp in vpathsums[up]:
                    sumval = vpathsums[up][vp]
                    if sumval == 0:
                        continue
                    for v2, vdiff, s in vpathdicts[index][vp]:
                        addsumval = s * sumval
                        for up2, udiff in newperms:
                            if vdiff + udiff == th[index]:
                                if up2 not in newpathsums:
                                    newpathsums[up2] = {}
                                newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + addsumval
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schub_coprod(perm, indices):
    mperm = perm
    indices = sorted(indices)
    ret_dict = {}
    k = len(indices)
    n = len(mperm)
    kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    max_required = max([kcd[i] + i for i in range(len(kcd))])
    kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
    N = len(kcd)
    kperm = ~(uncode(kcd2))
    coeff_dict = {kperm: 1}
    coeff_dict = schubmult(coeff_dict, mperm)

    inv_kperm = inv(kperm)
    inverse_kperm = ~kperm
    # total_sum = 0
    for perm, val in coeff_dict.items():
        if val == 0:
            continue
        # pperm = [*perm]
        downperm = perm*inverse_kperm
        if inv(downperm) == inv(perm) - inv_kperm:
            flag = True
            for i in range(N):
                if downperm[i] > N:
                    flag = False
                    break
            if not flag:
                continue
            firstperm = permtrim(list(downperm[0:N]))
            secondperm = permtrim([downperm[i] - N for i in range(N, len(downperm))])
            ret_dict[(firstperm, secondperm)] = val
    return ret_dict

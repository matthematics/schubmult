# from ._vars import (
#     var_y,
#     var_x,
#     var2,
#     var3,
#     q_var2,
# )
from functools import cached_property

from symengine import Add, Mul, Pow, expand, symarray

import schubmult.schubmult_double as norm_yz
from schubmult.perm_lib import (
    add_perm_dict,
    call_zvars,
    compute_vpathdicts,
    double_elem_sym_q,
    elem_sym_func_q,
    elem_sym_perms_q,
    elem_sym_perms_q_op,
    elem_sym_poly_q,
    inv,
    inverse,
    medium_theta,
    mulperm,
    permtrim,
    strict_theta,
    uncode,
)


class _gvars:
    @cached_property
    def n(self):
        return 100

    # @cached_property
    # def fvar(self):
    #     return 100

    @cached_property
    def var1(self):
        return tuple(symarray("x", self.n).tolist())

    @cached_property
    def var2(self):
        return tuple(symarray("y", self.n).tolist())

    @cached_property
    def var3(self):
        return tuple(symarray("z", self.n).tolist())

    @cached_property
    def q_var(self):
        return tuple(symarray("q", self.n).tolist())

    @cached_property
    def var_r(self):
        return symarray("r", 100)


_vars = _gvars()


def E(p, k, varl=_vars.var2[1:], var_x=_vars.var1):
    return elem_sym_poly_q(p, k, var_x[1:], varl)


def single_variable(coeff_dict, varnum, var2=_vars.var2, q_var=_vars.q_var):
    ret = {}
    for u in coeff_dict:
        if varnum - 1 < len(u):
            ret[u] = ret.get(u, 0) + var2[u[varnum - 1]] * coeff_dict[u]
        else:
            ret[u] = ret.get(u, 0) + var2[varnum] * coeff_dict[u]
        new_perms_k = elem_sym_perms_q(u, 1, varnum, q_var)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = elem_sym_perms_q(u, 1, varnum - 1, q_var)
        for perm, udiff, mul_val in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u] * mul_val
        for perm, udiff, mul_val in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u] * mul_val
    return ret


def mult_poly(coeff_dict, poly, var_x=_vars.var1, var_y=_vars.var2, q_var=_vars.q_var):
    if poly in var_x:
        return single_variable(coeff_dict, var_x.index(poly), var_y, q_var)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly(ret, a, var_x, var_y, q_var)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly(ret, base, var_x, var_y, q_var)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly(coeff_dict, a, var_x, var_y, q_var))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def nil_hecke(perm_dict, v, n, var2=_vars.var2, var3=_vars.var3):
    if v == (1, 2):
        return perm_dict
    th = strict_theta(inverse(v))
    mu = permtrim(uncode(th))
    vmu = permtrim(mulperm([*v], mu))

    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu, True)
    for u, val in perm_dict.items():
        vpathsums = {u: {(1, 2): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                newperms = elem_sym_perms_q_op(up, mx_th, th[index], n)
                for up2, udiff, mul_val in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0) * mul_val
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                0,
                            ) + s * sumval * elem_sym_func_q(
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
        toget = tuple(vmu)
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def elem_sym_func_q_q(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2, q_var=_vars.q_var):
    newk = k - udiff
    if newk < vdiff:
        return 0
    if newk == vdiff:
        return 1
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
    return elem_sym_poly_q(newk - vdiff, newk, yvars, zvars, q_var)


def schubpoly_quantum(v, var_x=_vars.var1, var_y=_vars.var2, q_var=_vars.q_var, coeff=1):
    th = strict_theta(inverse(v))
    mu = permtrim(uncode(th))
    vmu = permtrim(mulperm([*v], mu))
    if len(th) == 0:
        return coeff
    while th[-1] == 0:
        th.pop()
    vpathdicts = compute_vpathdicts(th, vmu)
    vpathsums = {(1, 2): {(1, 2): coeff}}
    inv_mu = inv(mu)
    inv_vmu = inv(vmu)
    inv_u = 0
    ret_dict = {}
    for index in range(len(th)):
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
                q_var,
            )
            for up2, udiff, mul_val in newperms:
                if up2 not in newpathsums:
                    newpathsums[up2] = {}
                for v in vpathdicts[index]:
                    sumval = vpathsums[up].get(v, 0) * mul_val
                    if sumval == 0:
                        continue
                    for v2, vdiff, s in vpathdicts[index][v]:
                        newpathsums[up2][v2] = newpathsums[up2].get(
                            v2,
                            0,
                        ) + s * sumval * elem_sym_func_q_q(
                            th[index],
                            index + 1,
                            up,
                            up2,
                            v,
                            v2,
                            udiff,
                            vdiff,
                            var_x,
                            var_y,
                            q_var,
                        )
        vpathsums = newpathsums
    toget = tuple(vmu)
    ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict[(1, 2)]


def schubmult(perm_dict, v, var2=_vars.var2, var3=_vars.var3, q_var=_vars.q_var):
    if v == (1, 2):
        return perm_dict
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
                    q_var,
                )
                for up2, udiff, mul_val in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0) * mul_val
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                0,
                            ) + s * sumval * elem_sym_func_q(
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
        toget = tuple(vmu)
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schubmult_db(perm_dict, v, var2=_vars.var2, var3=_vars.var3, q_var=_vars.q_var):
    if v == (1, 2):
        return perm_dict
    th = medium_theta(inverse(v))
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    mu = permtrim(uncode(th))
    vmu = permtrim(mulperm([*v], mu))
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}

    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu, True)
    for u, val in perm_dict.items():
        inv_u = inv(u)
        vpathsums = {u: {(1, 2): val}}
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
                                esim1 = (
                                    elem_sym_func_q(
                                        th[index],
                                        index + 1,
                                        up,
                                        up1,
                                        v,
                                        v2,
                                        udiff1,
                                        vdiff2,
                                        var2,
                                        var3,
                                    )
                                    * mul_val1
                                    * s2
                                )
                                mulfac = sumval * esim1
                                if (up1, udiff1, mul_val1) not in newpathsums0:
                                    newpathsums0[(up1, udiff1, mul_val1)] = {}
                                # newpathsums0[(up1, udiff1, mul_val1
                                newpathsums0[(up1, udiff1, mul_val1)][v2] = newpathsums0[(up1, udiff1, mul_val1)].get(v2, 0) + mulfac

                    for up1, udiff1, mul_val1 in newpathsums0:
                        for v in vpathdicts[index + 1]:
                            sumval = newpathsums0[(up1, udiff1, mul_val1)].get(v, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff2, s2 in vpathdicts[index + 1][v]:
                                for up2, udiff2, mul_val2 in newperms[(up1, udiff1, mul_val1)]:
                                    esim1 = (
                                        elem_sym_func_q(
                                            th[index + 1],
                                            index + 2,
                                            up1,
                                            up2,
                                            v,
                                            v2,
                                            udiff2,
                                            vdiff2,
                                            var2,
                                            var3,
                                        )
                                        * mul_val2
                                        * s2
                                    )
                                    mulfac = sumval * esim1
                                    if up2 not in newpathsums:
                                        newpathsums[up2] = {}
                                    newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + mulfac
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
                            sumval = vpathsums[up].get(v, 0) * mul_val
                            if sumval == 0:
                                continue
                            for v2, vdiff, s in vpathdicts[index][v]:
                                newpathsums[up2][v2] = newpathsums[up2].get(
                                    v2,
                                    0,
                                ) + s * sumval * elem_sym_func_q(
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
        toget = tuple(vmu)
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def div_diff(v, w, var2=_vars.var2, var3=_vars.var3):
    coeff_dict = {v: 1}
    coeff_dict = norm_yz.schubmult_down(coeff_dict, w, var2, var3)
    return coeff_dict.get((1, 2), 0)


def sum_q_dict(q_dict1, q_dict2):
    ret = {**q_dict1}
    for key in q_dict2:
        ret[key] = ret.get(key, 0) + q_dict2[key]
    return ret


def mul_q_dict(q_dict1, q_dict2):
    ret = {}
    for key1 in q_dict1:
        for key2 in q_dict2:
            key3 = key1 * key2
            ret[key3] = ret.get(key3, 0) + q_dict1[key1] * q_dict2[key2]
    return ret


def factor_out_q_keep_factored(poly):
    ret = {}
    if str(poly).find("q") == -1:
        ret[1] = poly
        return ret
    if poly in _vars.q_var:
        ret[poly] = 1
        return ret
    if isinstance(poly, Add):
        ag = poly.args
        ret = factor_out_q_keep_factored(ag[0])
        for i in range(1, len(ag)):
            ret = sum_q_dict(ret, factor_out_q_keep_factored(ag[i]))
        return ret
    if isinstance(poly, Mul):
        ag = poly.args
        ret = factor_out_q_keep_factored(ag[0])
        for i in range(1, len(ag)):
            ret = mul_q_dict(ret, factor_out_q_keep_factored(ag[i]))
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])

        ret = factor_out_q_keep_factored(base)
        ret0 = dict(ret)
        for _ in range(exponent - 1):
            ret = mul_q_dict(ret, ret0)

        # print(f"exponent {exponent}")
        # work_val = factor_out_q_keep_factored(base)
        # ret = {1: 1}
        # while exponent > 0:
        # if exponent % 2 == 1:
        # if ret == {1: 1}:
        # ret = {**work_val}
        # else:
        # ret = mul_q_dict(ret,work_val)
        # exponent -= 1
        # else:
        # work_val = mul_q_dict(work_val,work_val)
        # exponent //= 2
        return ret
    return ret


def factor_out_q(poly):
    coeff_dict = expand(poly).as_coefficients_dict()
    ret = {}
    for key in coeff_dict:
        coeff = coeff_dict[key]
        if coeff == 0:
            continue
        q_part = 1
        yz_part = coeff
        if isinstance(key, Mul):
            for var_maybe_pow in key.args:
                if isinstance(var_maybe_pow, Pow):
                    real_var = var_maybe_pow.args[0]
                    if real_var in _vars.q_var:
                        q_part *= var_maybe_pow
                    else:
                        yz_part *= var_maybe_pow
                else:
                    real_var = var_maybe_pow
                    if real_var in _vars.q_var:
                        q_part *= var_maybe_pow
                    else:
                        yz_part *= var_maybe_pow
        elif isinstance(key, Pow):
            real_var = key.args[0]
            if real_var in _vars.q_var:
                q_part *= key
            else:
                yz_part *= key
        elif key in _vars.q_var:
            q_part *= key
        else:
            yz_part *= key

        ret[q_part] = ret.get(q_part, 0) + yz_part
    return ret

# from ._vars import (
#     var_y,
#     var_x,
#     var2,
#     var3,
#     q_var2,
# )
from functools import cache

import numpy as np

import schubmult.mult.double as norm_yz
import schubmult.mult.positivity as pos
from schubmult.rings.poly_lib import _vars, call_zvars, elem_sym_func_q, elem_sym_poly_q, q_vector
from schubmult.rings.variables import CustomGeneratingSet, GeneratingSet_base
from schubmult.schub_lib.perm_lib import Permutation, code, inv, longest_element, medium_theta, strict_theta, uncode
from schubmult.symbolic import Add, Mul, Pow, S, expand, sympify
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import (
    add_perm_dict,
    count_less_than,
    is_parabolic,
    omega,
)
from schubmult.utils.schub_lib import check_blocks, compute_vpathdicts, double_elem_sym_q, elem_sym_perms_q, elem_sym_perms_q_op, elem_sym_positional_perms_q, reduce_q_coeff

logger = get_logger(__name__)


# def E(p, k, varl=None, var_x=None):
#     return elem_sym_poly_q(p, k, var_x[1:], varl)


def single_variable(coeff_dict, varnum, var_y=None, q_var=_vars.q_var):
    ret = {}
    for u in coeff_dict:
        ret[u] = ret.get(u, S.Zero) + var_y[u[varnum - 1]] * coeff_dict[u]
        # else:
        #     ret[u] = ret.get(u, 0) + var_y[varnum] * coeff_dict[u]
        new_perms_k = elem_sym_positional_perms_q(u, 1, varnum, q_var=q_var)
        # new_perms_km1 = []
        # if varnum > 1:
        #     new_perms_km1 = elem_sym_perms_q(u, 1, varnum - 1, q_var)
        for perm, udiff, sign, mul_val in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, S.Zero) + coeff_dict[u] * mul_val * sign
        # for perm, udiff, mul_val in new_perms_km1:
        #     if udiff == 1:
        #         ret[perm] = ret.get(perm, 0) - coeff_dict[u] * mul_val
    return ret


def mult_poly_q_double(coeff_dict, poly, var_x=None, var_y=None, q_var=_vars.q_var):
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    # logger.debug(f"{poly=} {list(var_x)=}")
    if var_x.index(poly) != -1:
        # logger.debug(f"yay {var_x.index(poly)=}")
        return single_variable(coeff_dict, var_x.index(poly), var_y, q_var)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_q_double(ret, a, var_x, var_y, q_var)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_q_double(ret, base, var_x, var_y, q_var)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_q_double(coeff_dict, a, var_x, var_y, q_var))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def mult_poly_q_double_alt(coeff_dict, poly, var_x=None, var_y=None, q_var=_vars.q_var):
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    if var_x.index(poly) != -1:
        return single_variable(coeff_dict, var_x.index(poly), var_y)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            s_d = mult_poly_q_double_alt({Permutation([]): S.One}, a, var_x, var_y, q_var)
            ret = schubmult_q_double_dict_fast(ret, s_d, var_y, var_y, q_var)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        s_d = mult_poly_q_double_alt({Permutation([]): S.One}, base, var_x, var_y, q_var)
        for i in range(int(exponent)):
            ret = schubmult_q_double_dict_fast(ret, s_d, var_y, var_y, q_var)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_q_double_alt(coeff_dict, a, var_x, var_y, q_var))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def nil_hecke(perm_dict, v, n, var2=None, var3=None):
    if v == Permutation([1, 2]):
        return perm_dict
    th = strict_theta(~v)
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
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


@cache
def schubmult_q_double_pair(perm1, perm2, var2=None, var3=None, q_var=None):
    return schubmult_q_double_fast({perm1: 1}, perm2, var2, var3, q_var)


@cache
def schubmult_q_double_pair_generic(perm1, perm2):
    return schubmult_q_double_fast({perm1: 1}, perm2, _vars.var_g1, _vars.var_g2, _vars.q_var)


@cache
def schubmult_q_generic_partial_posify(u2, v2):
    # logger.debug("Line number")
    return {w2: q_partial_posify_generic(val, u2, v2, w2) for w2, val in schubmult_q_double_pair_generic(u2, v2).items()}


def q_posify(u, v, w, val, var2, var3, q_var, msg):
    # logger.debug(f"Line number {val=} {u=} {v=} {w=}")
    if not v.has_pattern([1, 4, 3, 2]) and not v.has_pattern([3, 1, 2]):
        return schubmult_q_double_fast(u, v, var2, var3, q_var).get(w, S.Zero)
    try:
        val2 = int(expand(val))
    except Exception:
        # logger.debug("Line number")
        val2 = 0
        q_dict = factor_out_q_keep_factored(val)
        # logger.debug(f"{q_dict=}")
        # logger.debug("Line number")
        for q_part in q_dict:
            try:
                val2 += q_part * int(q_dict[q_part])
            except Exception:
                try:
                    # logger.debug("Line number")
                    if code(~v) == medium_theta(~v):
                        val2 += q_part * q_dict[q_part]
                    else:
                        q_part2 = q_part
                        qv = q_vector(q_part)
                        u2, v2, w2 = u, v, w
                        u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                        while did_one:
                            u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                        q_part2 = np.prod(
                            [q_var[i + 1] ** qv[i] for i in range(len(qv))],
                        )
                        if q_part2 == 1:
                            # reduced to classical coefficient
                            # logger.debug(f"{u=} {v=} {w=} {u2=} {v2=} {w2=} {q_part=} {q_dict[q_part]=}")
                            val2 += q_part * pos.posify(
                                q_dict[q_part],
                                u2,
                                v2,
                                w2,
                                var2,
                                var3,
                                msg,
                                False,
                            )
                        else:
                            val2 += q_part * pos.compute_positive_rep(
                                q_dict[q_part],
                                var2,
                                var3,
                                msg,
                            )
                            if val2 is None:
                                raise Exception
                except Exception:
                    import traceback

                    traceback.print_exc()
        if expand(val - val2) != 0:
            # logger.debug("Different")
            raise Exception
    return val2


# def q_posify(u, v, w, val, var2, var3, q_var, msg):
#     if expand(val) != 0:
#         try:
#             int(val)
#         except Exception:
#             val2 = 0
#             q_dict = factor_out_q_keep_factored(val)
#             for q_part in q_dict:
#                 try:
#                     val2 += q_part * int(q_dict[q_part])
#                 except Exception:
#                     # if same:
#                     #     to_add = q_part * expand(sympify(q_dict[q_part]).xreplace(subs_dict2))
#                     #     val2 += to_add
#                     # else:
#                     try:
#                         if code(~v) == medium_theta(~v):
#                             val2 += q_part * q_dict[q_part]
#                         else:
#                             q_part2 = q_part
#                             qv = q_vector(q_part)
#                             u2, v2, w2 = u, v, w
#                             u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
#                             while did_one:
#                                 u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
#                             q_part2 = np.prod(
#                                 [q_var[i + 1] ** qv[i] for i in range(len(qv))],
#                             )
#                             if q_part2 == 1:
#                                 # reduced to classical coefficient
#                                 val2 += q_part * pos.posify(
#                                     q_dict[q_part],
#                                     u2,
#                                     v2,
#                                     w2,
#                                     var2,
#                                     var3,
#                                     msg,
#                                     False,
#                                 )
#                             else:
#                                 val2 += q_part * pos.compute_positive_rep(
#                                     q_dict[q_part],
#                                     var2,
#                                     var3,
#                                     msg,
#                                     False,
#                                 )
#                     except Exception as e:
#                         # print(f"Exception: {e}")
#                         import traceback

#                         traceback.print_exc()
#                         exit(1)
#             if expand(val - val2) != 0:
#                raise Exception
#             val = val2
#         return val
#     return 0


def old_q_posify(u, v, w, val, var2, var3, q_var, msg):
    val2 = 0
    q_dict = factor_out_q_keep_factored(val)
    for q_part in q_dict:
        try:
            val2 += q_part * int(q_dict[q_part])
        except Exception:
            try:
                q_part2 = q_part
                qv = q_vector(q_part)
                u2, v2, w2 = u, v, w
                u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                while did_one:
                    u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                q_part2 = np.prod(
                    [q_var[i + 1] ** qv[i] for i in range(len(qv))],
                )
                if q_part2 == 1:
                    # reduced to classical coefficient
                    val2 += q_part * pos.posify(
                        q_dict[q_part],
                        u2,
                        v2,
                        w2,
                        var2,
                        var3,
                        msg,
                        False,
                    )
                else:
                    val2 += q_part * pos.compute_positive_rep(
                        q_dict[q_part],
                        var2,
                        var3,
                        msg,
                        False,
                    )
            except Exception:
                # print(f"Exception: {e}")
                import traceback

                traceback.print_exc()
                exit(1)
    if expand(val - val2) != 0:
        raise Exception
    return val2


def q_partial_posify_generic(val, u, v, w):
    if not v.has_pattern([1, 4, 3, 2]) and not v.has_pattern([3, 1, 2]):
        return val
    try:
        val2 = int(expand(val))
    except Exception:
        val2 = 0
        # logger.debug(f"{val=}")
        q_dict = factor_out_q_keep_factored(val)
        # logger.debug(f"{q_dict=}")
        for q_part in q_dict:
            try:
                val2 += q_part * int(q_dict[q_part])
            except Exception:
                try:
                    if code(~v) == medium_theta(~v):
                        val2 += q_part * q_dict[q_part]
                    else:
                        q_part2 = q_part
                        qv = q_vector(q_part)
                        u2, v2, w2 = u, v, w
                        u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                        while did_one:
                            u2, v2, w2, qv, did_one = reduce_q_coeff(u2, v2, w2, qv)
                        q_part2 = np.prod(
                            [_vars.q_var[i + 1] ** qv[i] for i in range(len(qv))],
                        )
                        if q_part2 == 1:
                            # reduced to classical coefficient
                            # logger.debug(f"{u=} {v=} {w=} {u2=} {v2=} {w2=} {q_part=} {q_dict[q_part]=}")
                            val2 += q_part * pos.posify_generic_partial(
                                q_dict[q_part],
                                u2,
                                v2,
                                w2,
                            )
                        else:
                            val2 += q_part * q_dict[q_part]
                except Exception as e:
                    logger.debug(f"Exception: {e}")

                    # import traceback

                    # traceback.print_exc()
        if expand(val - val2) != 0:
            raise Exception
    return val2


def apply_peterson_woodward(coeff_dict, parabolic_index, q_var=_vars.q_var):
    max_len = parabolic_index[-1] + 1
    w_P = longest_element(parabolic_index)
    w_P_prime = Permutation([1, 2])
    coeff_dict_update = {}
    for w_1 in coeff_dict.keys():
        val = coeff_dict[w_1]
        q_dict = factor_out_q_keep_factored(val)
        for q_part in q_dict:
            qv = q_vector(q_part)
            w = w_1
            good = True
            parabolic_index2 = []
            for i in range(len(parabolic_index)):
                if omega(parabolic_index[i], qv) == 0:
                    parabolic_index2 += [parabolic_index[i]]
                elif omega(parabolic_index[i], qv) != -1:
                    good = False
                    break
            if not good:
                continue
            w_P_prime = longest_element(parabolic_index2)
            if not check_blocks(qv, parabolic_index):
                continue
            w = (w * w_P_prime) * w_P
            if not is_parabolic(w, parabolic_index):
                continue

            if len(w) > max_len:
                continue
            new_q_part = np.prod(
                [q_var[index + 1 - count_less_than(parabolic_index, index + 1)] ** qv[index] for index in range(len(qv)) if index + 1 not in parabolic_index],
            )
            try:
                new_q_part = int(new_q_part)
            except Exception:
                pass
            q_val_part = q_dict[q_part]
            coeff_dict_update[w] = coeff_dict_update.get(w, 0) + new_q_part * q_val_part
    return coeff_dict_update


def elem_sym_func_q_q(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2, q_var=_vars.q_var):
    newk = k - udiff
    if newk < vdiff:
        return 0
    if newk == vdiff:
        return 1
    yvars = []
    # mlen = max(len(u1), len(u2))
    # u1 = [*u1] + [a + 1 for a in range(len(u1), mlen)]
    # u2 = [*u2] + [a + 1 for a in range(len(u2), mlen)]
    for j in range(k):
        if u1[j] == u2[j]:
            yvars += [varl1[u2[j]]]
    # for j in range(len(u1), min(k, len(u2))):
    #     if u2[j] == j + 1:
    #         yvars += [varl1[u2[j]]]
    # for j in range(len(u2), k):
    #     yvars += [varl1[j + 1]]
    zvars = [varl2[a] for a in call_zvars(v1, v2, k, i)]
    return elem_sym_poly_q(newk - vdiff, newk, yvars, zvars, q_var)


def schubpoly_quantum(v, var_x=None, var_y=None, q_var=_vars.q_var, coeff=1):
    th = strict_theta(~v)
    mu = uncode(th)
    vmu = v * mu
    if len(th) == 0:
        return coeff
    while len(th) > 0 and th[-1] == 0:
        th.pop()
    vpathdicts = compute_vpathdicts(th, vmu)
    vpathsums = {Permutation([1, 2]): {Permutation([1, 2]): coeff}}
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
    toget = vmu
    ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict[Permutation([1, 2])]


def schubmult_q_double(perm_dict, v, var2=None, var3=None, q_var=_vars.q_var):
    if v == Permutation([1, 2]):
        return perm_dict
    th = strict_theta(~v)
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}
    if len(th) == 0:
        return perm_dict
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
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schubmult_q_double_dict_fast(perm_dict1, perm_dict2, var2=None, var3=None, q_var=_vars.q_var):
    ret = {}
    for k, v in perm_dict2.items():
        ret = add_perm_dict(ret, {k2: v2 * v for k2, v2 in schubmult_q_double_fast(perm_dict1, k, var2, var3, q_var).items()})
    return ret


def schubmult_q_double_fast(perm_dict, v, var2=None, var3=None, q_var=_vars.q_var):
    if v == Permutation([1, 2]):
        return perm_dict
    th = medium_theta(~v)
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = inv(vmu)
    inv_mu = inv(mu)
    ret_dict = {}

    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        inv_u = inv(u)
        vpathsums = {u: {Permutation([]): val}}
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
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict



def div_diff(v, w, var2=None, var3=None):
    coeff_dict = {v: 1}
    coeff_dict = norm_yz.schubmult_down(coeff_dict, w, var2, var3)
    return coeff_dict.get(Permutation([1, 2]), 0)


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


def factor_out_q_keep_factored(poly, q_var=_vars.q_var):
    ret = {}
    # if str(poly).find("q") == -1:
    #     ret[1] = poly
    #     return ret
    if not isinstance(q_var, GeneratingSet_base):
        q_var = CustomGeneratingSet(q_var)
    # logger.debug(f"{poly=}")
    found_one = False
    for s in sympify(poly).free_symbols:
        if q_var.index(s) != -1:
            found_one = True
            # logger.debug("frobble bagel")

    if not found_one:
        ret[1] = poly
        return ret
    if q_var.index(poly) != -1:
        # logger.debug("it might be poke")
        # logger.debug(f"{poly=}")
        ret[poly] = 1
        return ret
    if isinstance(poly, Add):
        ag = poly.args
        ret = factor_out_q_keep_factored(ag[0], q_var)
        for i in range(1, len(ag)):
            ret = sum_q_dict(ret, factor_out_q_keep_factored(ag[i], q_var))
        return ret
    if isinstance(poly, Mul):
        ag = poly.args
        ret = factor_out_q_keep_factored(ag[0], q_var)
        for i in range(1, len(ag)):
            ret = mul_q_dict(ret, factor_out_q_keep_factored(ag[i], q_var))
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])

        ret = factor_out_q_keep_factored(base, q_var)
        ret0 = dict(ret)
        for _ in range(exponent - 1):
            ret = mul_q_dict(ret, ret0)
        return ret
    raise ValueError


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

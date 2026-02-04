from functools import cached_property

import schubmult.rings.variables as spl
import schubmult.utils.schub_lib as sss
from schubmult.schub_lib.permutation import (
    Permutation,
    uncode,
)
from schubmult.symbolic import Add, Mul, Pow
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

logger = get_logger(__name__)


class _gvars:
    @cached_property
    def n(self):
        return 100

    @cached_property
    def var_x(self):
        return spl.GeneratingSet("x")

    @cached_property
    def q_var(self):
        return spl.GeneratingSet("q")


_vars = _gvars()


def single_variable(coeff_dict, varnum, var_q=_vars.q_var):
    ret = {}
    for u in coeff_dict:
        new_perms_k = sss.elem_sym_perms_q(u, 1, varnum, var_q)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = sss.elem_sym_perms_q(u, 1, varnum - 1, var_q)
        # print(f"{new_perms_k=}")
        for perm, udiff, mul_val in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u] * mul_val
        for perm, udiff, mul_val in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u] * mul_val
    return ret


def mult_poly_q(coeff_dict, poly, var_x=_vars.var_x, var_q=_vars.q_var):
    if not isinstance(var_x, spl.GeneratingSet_base):
        var_x = spl.CustomGeneratingSet(var_x)
    # logger.debug(f"{poly=} {type(poly)=} {list(var_x)}")
    # logger.debug(f"{[type(v) for v in var_x]}")
    if var_x.index(poly) != -1:
        # logger.debug(f"Found {var_x.index(poly)=}")
        # print("bang")
        return single_variable(coeff_dict, var_x.index(poly), var_q=var_q)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_q(ret, a, var_x, var_q=var_q)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_q(ret, base, var_x, var_q=var_q)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_q(coeff_dict, a, var_x, var_q=var_q))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def schubmult_q_fast(perm_dict, v, q_var=_vars.q_var):
    if v.inv == 0:
        return perm_dict
    th = (~v).medium_theta()
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}

    thL = len(th)
    # if thL!=2 and len(set(thL))!=1:
    # raise ValueError("Not what I can do")
    vpathdicts = sss.compute_vpathdicts(th, vmu)
    # print(f"{vpathdicts=}")
    for u, val in perm_dict.items():
        inv_u = u.inv
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
                    inv_up = up.inv
                    newperms = sss.double_elem_sym_q(up, mx_th, mx_th1, th[index], q_var)
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
                    inv_up = up.inv
                    newperms = sss.elem_sym_perms_q(
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


def schubmult_q(perm_dict, v):
    th = (~v).strict_theta()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = sss.compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {u: {Permutation([]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = sss.elem_sym_perms_q(
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
        toget = vmu
        ret_dict = add_perm_dict(
            {ep: vpathsums[ep].get(toget, 0) for ep in vpathsums},
            ret_dict,
        )
    return ret_dict

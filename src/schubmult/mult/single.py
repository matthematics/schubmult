from functools import cached_property

from schubmult.rings.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
from schubmult.schub_lib.perm_lib import (
    Permutation,
    theta,
    uncode,
)
from schubmult.symbolic import Add, Mul, Pow
from schubmult.utils.logging import get_logger, init_logging
from schubmult.utils.perm_utils import (
    add_perm_dict,
)
from schubmult.utils.schub_lib import (
    compute_vpathdicts,
    elem_sym_perms,
    elem_sym_perms_op,
)

init_logging(debug=False)
logger = get_logger(__name__)

class _gvars:
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


# TODO: if need indexes, CustomGeneratingSet
def mult_poly_py(coeff_dict, poly, var_x=_vars.var_x):
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)

    if var_x.index(poly) != -1:
        # print(f"{poly=} {var_x._symbols_arr=} {var_x._symbols_arr.index(poly)=}")
        return single_variable(coeff_dict, var_x.index(poly))

    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_py(ret, a, var_x)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_py(ret, base, var_x)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_py(coeff_dict, a, var_x))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def schubmult_py(perm_dict, v):
    v = Permutation(v)
    # print(f"{v=}")
    vn1 = ~v
    th = theta(vn1)
    if len(th) == 0 or th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
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
        inv_u = u.inv
        vpathsums = {Permutation(u): {Permutation([1, 2]): val}}

        for index in range(len(th)):
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms(
                    up,
                    min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                    th[index],
                )
                # print(f"{up=}")
                for vp in vpathsums[up]:
                    # print(f"{vp=} {type(vp)=} {hash(vp)=}")
                    # print(f"{vpathsums[up]=} {vpathdicts[index]=}")
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

def schubmult_py_down(perm_dict, v):
    v = Permutation(v)
    # print(f"{v=}")
    vn1 = ~v
    th = theta(vn1)
    if len(th) == 0 or th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
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
        inv_u = u.inv
        vpathsums = {Permutation(u): {Permutation([1, 2]): val}}

        for index in range(len(th)):
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms_op(
                    up,
                    min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                    th[index],
                )
                # print(f"{up=}")
                for vp in vpathsums[up]:
                    # print(f"{vp=} {type(vp)=} {hash(vp)=}")
                    # print(f"{vpathsums[up]=} {vpathdicts[index]=}")
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


def schub_coprod_py(perm, indices):
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
    #logger.debug(f"{kperm.code}*{mperm.code}")
    coeff_dict = schubmult_py(coeff_dict, mperm)

    inv_kperm = kperm.inv
    inverse_kperm = ~kperm
    # total_sum = 0
    for perm, val in coeff_dict.items():
        if val == 0:
            continue
        # pperm = [*perm]
        downperm = perm * inverse_kperm
        if downperm.inv == perm.inv - inv_kperm:
            flag = True
            for i in range(N):
                if downperm[i] > N:
                    flag = False
                    break
            if not flag:
                continue
            firstperm = downperm[0:N]
            secondperm = Permutation([downperm[i] - N for i in range(N, len(downperm))])
            ret_dict[(firstperm, secondperm)] = val
    return ret_dict

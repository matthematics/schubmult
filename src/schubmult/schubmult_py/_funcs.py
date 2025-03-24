from functools import cached_property

from symengine import Add, Mul, Pow, symarray

from schubmult.perm_lib import (
    add_perm_dict,
    compute_vpathdicts,
    elem_sym_perms,
    inv,
    inverse,
    mulperm,
    permtrim,
    theta,
    uncode,
)


class _gvars:
    @cached_property
    def n(self):
        return 100

    @cached_property
    def var_x(self):
        return tuple(symarray("x", self.n).tolist())


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
    if poly in var_x:
        return single_variable(coeff_dict, var_x.index(poly))
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
    vn1 = inverse(v)
    th = theta(vn1)
    if th[0] == 0:
        return perm_dict
    mu = permtrim(uncode(th))
    vmu = permtrim(mulperm(list(v), mu))
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
        vpathsums = {u: {(1, 2): val}}

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
        toget = tuple(vmu)
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schub_coprod(perm, indices):
    mperm = tuple(perm)
    indices = sorted(indices)
    ret_dict = {}
    k = len(indices)
    n = len(mperm)
    kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    max_required = max([kcd[i] + i for i in range(len(kcd))])
    kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
    N = len(kcd)
    kperm = permtrim(inverse(uncode(kcd2)))
    coeff_dict = {tuple(kperm): 1}
    coeff_dict = schubmult(coeff_dict, tuple(mperm))

    inv_kperm = inv(kperm)
    inverse_kperm = inverse(kperm)
    # total_sum = 0
    for perm, val in coeff_dict.items():
        if val == 0:
            continue
        pperm = [*perm]
        downperm = mulperm(pperm, inverse_kperm)
        if inv(downperm) == inv(pperm) - inv_kperm:
            flag = True
            for i in range(N):
                if downperm[i] > N:
                    flag = False
                    break
            if not flag:
                continue
            firstperm = tuple(permtrim(list(downperm[0:N])))
            secondperm = tuple(permtrim([downperm[i] - N for i in range(N, len(downperm))]))
            ret_dict[(firstperm, secondperm)] = val
    return ret_dict

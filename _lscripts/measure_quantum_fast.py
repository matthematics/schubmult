"""Fast factor-only measurement of the quantum kernel (no symbolic verification).

Confirms at larger n:
  qmonos_max : max #q-monomials in any coefficient  (expect 1 => coeffs are monomials)
  front_max  : max #(up,v-label) entries per step
  maxlen_gap : max (ell(up) - (l_u+l_v))  (expect <= 0)
  Bsingle    : max single-column elem_sym_perms_q branch
  Bdouble    : max double_elem_sym_q total covers
  L_eff_ok   : whether L_eff <= min(n-1, l_v) always
  qdeg_max   : max q-degree
"""

import itertools
import sys
from math import comb

from schubmult import Permutation
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S, expand, sympify
from schubmult.utils import schub_lib as sss


def qmono_count(coeff):
    if coeff == 0 or coeff == 1 or coeff == S.One:
        return 1 if coeff != 0 else 0
    e = expand(sympify(coeff))
    if e == 0:
        return 0
    return len(e.args) if e.is_Add else 1


def run(u, v, stats):
    perm_dict = {u: S.One}
    if v.inv == 0:
        return
    th = (~v).medium_theta()
    if len(th) == 0:
        return
    while th[-1] == 0:
        th.pop()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    thL = len(th)
    vpathdicts = sss.compute_vpathdicts(th, vmu)
    L_eff = 0
    for uu, val in perm_dict.items():
        inv_u = uu.inv
        vpathsums = {uu: {Permutation([1, 2]): val}}
        for index in range(thL):
            if index > 0 and th[index - 1] == th[index]:
                continue
            L_eff += 1
            fsize = sum(len(vpathsums[up]) for up in vpathsums)
            stats["front"] = max(stats["front"], fsize)
            for up in vpathsums:
                stats["lgap"] = max(stats["lgap"], up.inv - (u.inv + v.inv))
                for lab, cf in vpathsums[up].items():
                    stats["qmonos"] = max(stats["qmonos"], qmono_count(cf))
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
                    newperms = sss.double_elem_sym_q(up, mx_th, mx_th1, th[index])
                    stats["Bdouble"] = max(stats["Bdouble"], sum(len(newperms[k]) for k in newperms))
                    for v_ in vpathdicts[index]:
                        sumval = vpathsums[up].get(v_, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff2, s2 in vpathdicts[index][v_]:
                            for up1, udiff1, mul_val1 in newperms:
                                if (up1, udiff1, mul_val1) not in newpathsums0:
                                    newpathsums0[(up1, udiff1, mul_val1)] = {}
                                if udiff1 + vdiff2 == th[index]:
                                    newpathsums0[(up1, udiff1, mul_val1)][v2] = newpathsums0[(up1, udiff1, mul_val1)].get(v2, 0) + s2 * sumval * mul_val1
                    for up1, udiff1, mul_val1 in newpathsums0:
                        for v_ in vpathdicts[index + 1]:
                            sumval = newpathsums0[(up1, udiff1, mul_val1)].get(v_, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff2, s2 in vpathdicts[index + 1][v_]:
                                for up2, udiff2, mul_val2 in newperms[(up1, udiff1, mul_val1)]:
                                    if up2 not in newpathsums:
                                        newpathsums[up2] = {}
                                    if udiff2 + vdiff2 == th[index + 1]:
                                        newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s2 * sumval * mul_val2
            else:
                newpathsums = {}
                for up in vpathsums:
                    inv_up = up.inv
                    newperms = sss.elem_sym_perms_q(up, min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu), th[index])
                    stats["Bsingle"] = max(stats["Bsingle"], len(newperms))
                    for up2, udiff, mul_val in newperms:
                        if up2 not in newpathsums:
                            newpathsums[up2] = {}
                        for v_ in vpathdicts[index]:
                            sumval = vpathsums[up].get(v_, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff, s in vpathdicts[index][v_]:
                                if udiff + vdiff == th[index]:
                                    newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s * sumval * mul_val
            vpathsums = newpathsums
    lbound = min(len(u) - 1, v.inv)
    if L_eff > lbound:
        stats["Leff_viol"] = max(stats["Leff_viol"], (L_eff, lbound, list(u), list(v)))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for n in range(2, nmax + 1):
        stats = {"front": 0, "lgap": -999, "qmonos": 0, "Bsingle": 0, "Bdouble": 0, "Leff_viol": (0, 0, None, None)}
        for u, v in itertools.product(Permutation.all_permutations(n), repeat=2):
            run(u, v, stats)
        print(
            f"n={n}: front={stats['front']} lgap={stats['lgap']} qmonos={stats['qmonos']} "
            f"Bsingle={stats['Bsingle']} Bdouble={stats['Bdouble']} "
            f"(2^n-1)={2**n-1} (2^n-1)^2={(2**n-1)**2} Leff_viol={stats['Leff_viol'][:2]}",
            flush=True,
        )


if __name__ == "__main__":
    main()

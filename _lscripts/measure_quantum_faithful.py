"""Faithful instrumented replica of schubmult_q_fast.

Goal: verify the replica matches the real kernel EXACTLY over all (u,v) up to n,
then measure the true quantum complexity factors:

  L_eff     = number of processed layer-steps (a paired double-column counts once)
  front     = max over steps of #(up-perm, v-label) coefficient entries
  maxlen    = max length ell(up) of any front permutation reached
  Bsingle   = max single-column elem_sym_perms_q branch actually used
  Bdouble   = max double_elem_sym_q total covers actually used
  qmonos    = max #q-monomials in any single coefficient (q-polynomial size)
  qdeg      = max total q-degree appearing in any coefficient
  Delta     = kdown-triple (v-side) factor  [UNCHANGED from classical]
"""

import itertools
import sys

from schubmult import Permutation
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S, expand, sympify
from schubmult.utils import schub_lib as sss
from schubmult.utils.perm_utils import add_perm_dict

from schubmult.mult.quantum import schubmult_q_fast


def qmono_count(coeff):
    e = expand(sympify(coeff))
    if e == 0:
        return 0
    try:
        args = e.args
        if e.is_Add:
            return len(args)
        return 1
    except Exception:
        return 1


def qdegree(coeff):
    e = sympify(coeff)
    try:
        d = 0
        # total degree in q-variables of each monomial; take max
        e = expand(e)
        terms = e.args if e.is_Add else [e]
        for t in terms:
            td = 0
            for f, expo in t.as_powers_dict().items():
                if str(f).startswith("q"):
                    td += int(expo)
            d = max(d, td)
        return d
    except Exception:
        return 0


def instrumented(u, v):
    stats = {"L_eff": 0, "front": 0, "maxlen": 0, "Bsingle": 0, "Bdouble": 0, "qmonos": 0, "qdeg": 0, "Delta": 0}
    perm_dict = {u: S.One}
    if v.inv == 0:
        return dict(perm_dict), stats
    th = (~v).medium_theta()
    if len(th) == 0:
        return dict(perm_dict), stats
    while th[-1] == 0:
        th.pop()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    thL = len(th)
    vpathdicts = sss.compute_vpathdicts(th, vmu)
    for uu, val in perm_dict.items():
        inv_u = uu.inv
        vpathsums = {uu: {Permutation([1, 2]): val}}
        for index in range(thL):
            if index > 0 and th[index - 1] == th[index]:
                continue
            # front measurement at start of a processed step
            fsize = sum(len(vpathsums[up]) for up in vpathsums)
            stats["front"] = max(stats["front"], fsize)
            for up in vpathsums:
                stats["maxlen"] = max(stats["maxlen"], up.inv)
                for lab, cf in vpathsums[up].items():
                    stats["qmonos"] = max(stats["qmonos"], qmono_count(cf))
                    stats["qdeg"] = max(stats["qdeg"], qdegree(cf))
            for vp in vpathdicts[index]:
                stats["Delta"] = max(stats["Delta"], len(vpathdicts[index][vp]))
            stats["L_eff"] += 1
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            if index < len(th) - 1 and th[index] == th[index + 1]:
                mx_th1 = 0
                for vp in vpathdicts[index + 1]:
                    for v2, vdiff, s in vpathdicts[index + 1][vp]:
                        mx_th1 = max(mx_th1, th[index + 1] - vdiff)
                for vp in vpathdicts[index + 1]:
                    stats["Delta"] = max(stats["Delta"], len(vpathdicts[index + 1][vp]))
                newpathsums = {}
                for up in vpathsums:
                    newpathsums0 = {}
                    inv_up = up.inv
                    newperms = sss.double_elem_sym_q(up, mx_th, mx_th1, th[index])
                    stats["Bdouble"] = max(stats["Bdouble"], sum(len(newperms[k]) for k in newperms))
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff2, s2 in vpathdicts[index][v]:
                            for up1, udiff1, mul_val1 in newperms:
                                if (up1, udiff1, mul_val1) not in newpathsums0:
                                    newpathsums0[(up1, udiff1, mul_val1)] = {}
                                if udiff1 + vdiff2 == th[index]:
                                    newpathsums0[(up1, udiff1, mul_val1)][v2] = newpathsums0[(up1, udiff1, mul_val1)].get(v2, 0) + s2 * sumval * mul_val1
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
                    newperms = sss.elem_sym_perms_q(up, min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu), th[index])
                    stats["Bsingle"] = max(stats["Bsingle"], len(newperms))
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
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict, stats


def norm(d):
    out = {}
    for k, val in d.items():
        e = expand(sympify(val))
        if e != 0:
            out[k] = e
    return out


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    mism = 0
    checked = 0
    agg = {"L_eff": 0, "front": 0, "maxlen": 0, "Bsingle": 0, "Bdouble": 0, "qmonos": 0, "qdeg": 0, "Delta": 0}
    argmax = {k: None for k in agg}
    worst_over = (0, None)  # maxlen - (l_u + l_v)
    for n in range(2, nmax + 1):
        for u, v in itertools.product(Permutation.all_permutations(n), repeat=2):
            got, stats = instrumented(u, v)
            ref = schubmult_q_fast({u: S.One}, v)
            if norm(got) != norm(ref):
                mism += 1
                if mism <= 3:
                    print("MISMATCH", list(u), list(v))
            checked += 1
            for k in agg:
                if stats[k] > agg[k]:
                    agg[k] = stats[k]
                    argmax[k] = (n, list(u), list(v))
            over = stats["maxlen"] - (u.inv + v.inv)
            if over > worst_over[0]:
                worst_over = (over, (n, list(u), list(v)))
    print(f"checked {checked} pairs up to n={nmax}, mismatches={mism}", flush=True)
    for k in agg:
        print(f"  max {k} = {agg[k]}  at {argmax[k]}", flush=True)
    print(f"  worst (maxlen - (l_u+l_v)) = {worst_over[0]} at {worst_over[1]}", flush=True)


if __name__ == "__main__":
    main()

"""Exact operation count of the schubmult inner loop.

Replicates the kernel of schubmult_py (src/schubmult/mult/single.py, lines 99-131)
verbatim in its loop structure, and counts the number of executions of the
innermost body -- the `if vdiff + udiff == th[index]` test. That count IS the
work W of the algorithm (up to an O(n) per-op factor for permutation hashing).

We record, per layer i:
  N_i    = number of (up, vp) nodes in the front
  beta_i = max over up of |elem_sym_perms(up,...)| (PIERI u-branching)
  delta_i= max over vp of |vpathdicts[i][vp]|      (KDOWN v-triples)
and the exact inner-body count
  W_i    = sum_{up} |newperms(up)| * sum_{vp in labels(up)} |vpathdicts[i][vp]|.

Then W = sum_i W_i. We compare W against candidate closed-form bounds.
"""

import itertools

from schubmult import Permutation, Sx
from schubmult.combinatorics.permutation import uncode
from schubmult.utils.schub_lib import compute_vpathdicts, elem_sym_perms


def exact_work(u, v):
    """Return (W, Wfilt, L, maxfront_nodes, max_beta, max_delta) for S_u * S_v."""
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return {"W": 0, "Wfilt": 0, "L": 0, "maxN": 1, "maxbeta": 0, "maxdelta": 0, "M": len(Sx(u) * Sx(v))}
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    th = list(th)
    while th and th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
    mx_th = [0] * thL
    for index in range(thL):
        for vp in vpathdicts[index]:
            for _v2, vdiff, _s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    inv_u = u.inv
    vpathsums = {u: {Permutation([1, 2]): 1}}

    W = 0        # total innermost-body executions (the vdiff+udiff test)
    Wfilt = 0    # executions that pass the filter (edges actually created)
    maxN = 1
    maxbeta = 0
    maxdelta = 0
    for index in range(thL):
        newpathsums = {}
        N = sum(len(vpathsums[up]) for up in vpathsums)  # (up,vp) nodes
        maxN = max(maxN, N)
        for up in vpathsums:
            inv_up = up.inv
            newperms = elem_sym_perms(up, min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)), th[index])
            b = len(newperms)
            maxbeta = max(maxbeta, b)
            for vp in vpathsums[up]:
                triples = vpathdicts[index][vp]
                maxdelta = max(maxdelta, len(triples))
                for _v2, vdiff, _s in triples:
                    for _up2, udiff in newperms:
                        W += 1
                        if vdiff + udiff == th[index]:
                            Wfilt += 1
                            newpathsums.setdefault(_up2, {})
                            newpathsums[_up2][_v2] = newpathsums[_up2].get(_v2, 0) + 1
        vpathsums = newpathsums

    return {"W": W, "Wfilt": Wfilt, "L": thL, "maxN": maxN, "maxbeta": maxbeta, "maxdelta": maxdelta, "M": len(Sx(u) * Sx(v))}


def main():
    import sys

    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    worstW = (0, None, None)
    worst_ratio_edges = (0, None, None)  # W / Wfilt
    total_checked = 0
    for n in range(2, nmax + 1):
        perms = Permutation.all_permutations(n)
        for u, v in itertools.product(perms, repeat=2):
            r = exact_work(u, v)
            total_checked += 1
            if r["W"] > worstW[0]:
                worstW = (r["W"], list(u), list(v), r)
            if r["Wfilt"] > 0:
                ratio = r["W"] / r["Wfilt"]
                if ratio > worst_ratio_edges[0]:
                    worst_ratio_edges = (ratio, list(u), list(v), r)
    print(f"checked {total_checked} pairs up to n={nmax}")
    print(f"worst W: {worstW[0]}  u={worstW[1]} v={worstW[2]}")
    print(f"   detail: {worstW[3]}")
    print(f"worst W/Wfilt ratio: {worst_ratio_edges[0]:.2f}  u={worst_ratio_edges[1]} v={worst_ratio_edges[2]}")
    print(f"   detail: {worst_ratio_edges[3]}")


if __name__ == "__main__":
    main()

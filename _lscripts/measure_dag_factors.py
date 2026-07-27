"""Measure the two DAG factors of the transition algorithm separately.

For each (u, v) we count:
  - Mu_tilde : positive up-path leaves on the u-side (elem_sym_perms),
               summed multiplicatively across the theta columns. This is the
               cancellation-free factor.
  - Mv_tilde : signed down-path leaves on the v-side (compute_vpathdicts /
               kdown_perms). This is the signed factor.
  - edges    : total (up x down) transition edges actually traversed = the
               real work of schubmult_double.
  - M        : true number of nonzero terms in the output.

The point: test whether edges ~ Mu_tilde * Mv_tilde (multiplicative), whether
Mu_tilde tracks M (unavoidable output), and how big Mv_tilde grows (the
putative removable inefficiency).
"""

import itertools

from schubmult import Permutation, Sx
from schubmult.combinatorics.permutation import uncode
from schubmult.utils.schub_lib import compute_vpathdicts, elem_sym_perms


def count_factors(u, v):
    """Return (Mu_tilde, Mv_tilde, edges, M) for the product S_u * S_v."""
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return (1, 1, 0, 1)
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    while th and th[-1] == 0:
        th.pop()
    thL = len(th)

    vpathdicts = compute_vpathdicts(th, vmu)

    # v-side width per column and total signed down-edges
    Mv_tilde = 0
    for index in range(thL):
        for vp in vpathdicts[index]:
            Mv_tilde += len(vpathdicts[index][vp])

    mx_th = [0] * thL
    for index in range(thL):
        for vp in vpathdicts[index]:
            for _v2, vdiff, _s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    inv_u = u.inv
    # Track the u-front width layer by layer (positive factor).
    vpathsums = {u: {Permutation([1, 2]): 1}}
    Mu_tilde = 0
    edges = 0
    for index in range(thL):
        newpathsums = {}
        for up in vpathsums:
            inv_up = up.inv
            newperms = elem_sym_perms(
                up,
                min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                th[index],
            )
            Mu_tilde += len(newperms)
            for vp in vpathsums[up]:
                for _v2, vdiff, _s in vpathdicts[index][vp]:
                    for up2, udiff in newperms:
                        if vdiff + udiff == th[index]:
                            edges += 1
                            newpathsums.setdefault(up2, {})
                            newpathsums[up2][_v2] = 1
        vpathsums = newpathsums

    result = Sx(u) * Sx(v)
    M = len(result)
    return (Mu_tilde, Mv_tilde, edges, M)


def main():
    import sys

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    perms = Permutation.all_permutations(n)
    print(f"n={n}  ({len(perms)} perms)")
    print(f"{'u':<18}{'v':<18}{'Mu~':>8}{'Mv~':>8}{'edges':>10}{'M':>8}{'edges/M':>10}{'Mu~/M':>8}")
    worst_mv = (0, None, None)
    for u, v in itertools.product(perms, repeat=2):
        Mu, Mv, edges, M = count_factors(u, v)
        if edges == 0:
            continue
        if Mv > worst_mv[0]:
            worst_mv = (Mv, u, v)
        # only print a sample: nontrivial v
        if v.inv >= n - 1 and u.inv <= 2:
            print(f"{str(list(u)):<18}{str(list(v)):<18}{Mu:>8}{Mv:>8}{edges:>10}{M:>8}{edges / M:>10.2f}{Mu / M:>8.2f}")
    print(f"\nworst Mv~: {worst_mv[0]} at u={list(worst_mv[1])} v={list(worst_mv[2])}")


if __name__ == "__main__":
    main()

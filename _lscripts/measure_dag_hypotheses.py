"""Focused test of the two hypotheses:

  H1 (Mu unavoidable): for v dominant, the v-factor is a single chain
     (Mv_tilde collapses) and Mu_tilde == M exactly (cancellation-free,
     leaves ARE the output).

  H2 (Mv is u-independent overhead): Mv_tilde depends only on v, not u.
     Measure Mv_tilde(v) and how it grows with l(v).
"""

import itertools

from schubmult import Permutation, Sx
from schubmult.combinatorics.permutation import uncode
from schubmult.utils.schub_lib import compute_vpathdicts, elem_sym_perms


def factors(u, v):
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return (1, 0, 0, 1)
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    while th and th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)

    Mv = 0
    for index in range(thL):
        for vp in vpathdicts[index]:
            Mv += len(vpathdicts[index][vp])

    mx_th = [0] * thL
    for index in range(thL):
        for vp in vpathdicts[index]:
            for _v2, vdiff, _s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    inv_u = u.inv
    vpathsums = {u: {Permutation([1, 2]): 1}}
    Mu = 0
    edges = 0
    for index in range(thL):
        newpathsums = {}
        for up in vpathsums:
            inv_up = up.inv
            newperms = elem_sym_perms(up, min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)), th[index])
            Mu += len(newperms)
            for vp in vpathsums[up]:
                for _v2, vdiff, _s in vpathdicts[index][vp]:
                    for up2, udiff in newperms:
                        if vdiff + udiff == th[index]:
                            edges += 1
                            newpathsums.setdefault(up2, {})
                            newpathsums[up2][_v2] = 1
        vpathsums = newpathsums

    M = len(Sx(u) * Sx(v))
    return (Mu, Mv, edges, M)


def is_dominant(v):
    # dominant = code weakly decreasing = 132-avoiding; theta == code
    c = list(v.code)
    return all(c[i] >= c[i + 1] for i in range(len(c) - 1))


def main():
    import sys

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    perms = Permutation.all_permutations(n)

    # H1: dominant v
    print("=== H1: v dominant  =>  Mv~ collapses, Mu~ == M ? ===")
    h1_ok = True
    dom_vs = [v for v in perms if is_dominant(v) and v.inv > 0]
    for v in dom_vs:
        for u in perms:
            Mu, Mv, edges, M = factors(u, v)
            if edges == 0:
                continue
            if Mu != M:
                h1_ok = False
                print(f"  Mu!=M: u={list(u)} v={list(v)} Mu={Mu} M={M} Mv={Mv}")
    print(f"H1 (Mu~==M for all dominant v, all u): {'HOLDS' if h1_ok else 'FAILS'} at n={n}")
    # report Mv range for dominant v
    domMv = {tuple(v): factors(Permutation([]), v)[1] for v in dom_vs}
    print(f"  dominant-v Mv~ values: min={min(domMv.values())} max={max(domMv.values())}")

    # H2: Mv~ depends only on v
    print("\n=== H2: Mv~ depends only on v (not u) ? ===")
    h2_ok = True
    sample_vs = [v for v in perms if v.inv >= 1]
    for v in sample_vs:
        vals = set()
        for u in perms:
            _, Mv, edges, _ = factors(u, v)
            if edges > 0:
                vals.add(Mv)
        if len(vals) > 1:
            h2_ok = False
            print(f"  Mv~ varies with u for v={list(v)}: {sorted(vals)}")
    print(f"H2 (Mv~ is u-independent): {'HOLDS' if h2_ok else 'FAILS'} at n={n}")

    # H2b: growth of Mv~ with l(v)
    print("\n=== Mv~ growth with l(v) ===")
    by_len = {}
    for v in sample_vs:
        _, Mv, _, _ = factors(Permutation([]), v)
        by_len.setdefault(v.inv, []).append(Mv)
    for lv in sorted(by_len):
        arr = by_len[lv]
        print(f"  l(v)={lv:>2}: Mv~ max={max(arr):>4}  mean={sum(arr) / len(arr):>6.2f}  (count {len(arr)})")


if __name__ == "__main__":
    main()

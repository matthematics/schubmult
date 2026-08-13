#!/usr/bin/env python3
"""Growth of the u-free ratios controlling output-polynomiality in the average.

The chain that needs to close is

    Lambda(mdom v)  <=  N^up  <=  poly(n, Ntilde),      Ntilde >= mean_{v' in C} Lambda(v')

so the decisive u-free quantity is  R(C) = Lambda(b) / mean_{v' in C} Lambda(v').
We report its max over classes for each n, together with the exponent
log Lambda(b) / log(mean Lambda) that a pure power bound would need.
"""

import math
import sys
from collections import defaultdict

from schubmult.combinatorics.permutation import Permutation


def downset_sizes(n):
    perms = Permutation.all_permutations(n)
    index = {w: i for i, w in enumerate(perms)}
    down = {}
    for w in sorted(perms, key=lambda w: w.inv):
        bits = 1 << index[w]
        wi = ~w
        for d in wi.zero_indexed_descents():
            bits |= down[~(wi.swap(d, d + 1))]
        down[w] = bits
    return {w: bin(b).count("1") for w, b in down.items()}


def sweep(n):
    lam = downset_sizes(n)
    buckets = defaultdict(list)
    for v in Permutation.all_permutations(n):
        buckets[(v.mul_sortable(), v.mul_dominant())].append(v)

    best_shadow = None
    best_exp = None
    best_star = None
    best_lamratio = None
    for (a, b), members in buckets.items():
        size = len(members)
        la, lb = lam[a], lam[b]
        mean_lam = sum(lam[v] for v in members) / size
        shadow = lb / mean_lam
        expo = math.log(lb) / math.log(mean_lam) if mean_lam > 1 else 1.0
        star = size / la
        lamratio = lb / la
        rec = (a, b, size, la, lb, mean_lam, shadow, expo, star, lamratio)
        if best_shadow is None or shadow > best_shadow[6]:
            best_shadow = rec
        if best_exp is None or expo > best_exp[7]:
            best_exp = rec
        if best_star is None or star > best_star[8]:
            best_star = rec
        if best_lamratio is None or lamratio > best_lamratio[9]:
            best_lamratio = rec
    return len(buckets), best_shadow, best_exp, best_star, best_lamratio


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    print(f"{'n':>2} {'#cls':>6} {'maxR=Lb/meanL':>14} {'max exponent':>13} {'max|C|/La':>10} {'max Lb/La':>10}")
    for n in range(3, hi + 1):
        ncls, sh, ex, st, lr = sweep(n)
        print(f"{n:>2} {ncls:>6} {sh[6]:>14.4f} {ex[7]:>13.4f} {st[8]:>10.4f} {lr[9]:>10.4f}")
        print(f"     worst R at class {sh[0]} .. {sh[1]}: |C|={sh[2]} Lambda(a)={sh[3]} Lambda(b)={sh[4]} mean={sh[5]:.2f}")

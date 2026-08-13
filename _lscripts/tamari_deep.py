#!/usr/bin/env python3
"""Push the shadow-ratio sweep beyond n=8.

Two accelerations:

1. b = mdom(v) has b^{-1} dominant (= 132-avoiding), so Lambda(b) = |[e,b]| is
   given in closed form by Knuth's hook length formula for forests,
   Lambda(b) = N! / prod_i |down(i)|,  down(i) = {j <= i : sigma(j) <= sigma(i)},
   for sigma = b^{-1}.  This is O(n^2) per class instead of a downset traversal.

2. Every member of the class lies below b, so once Lambda(b) is small the whole
   class average can be computed by a bitset DP restricted to the ideal [e,b],
   never touching the rest of S_n.  Classes with Lambda(b) above a cap are
   skipped and reported separately.
"""

import math
import sys
from collections import defaultdict

from schubmult.combinatorics.permutation import Permutation


def hook_lambda(b):
    """|[e,b]| via the forest hook length formula, valid when b^{-1} is dominant."""
    sigma = ~b
    N = len(sigma)
    prod = 1
    for i in range(1, N + 1):
        prod *= 1 + sum(1 for j in range(1, i) if sigma[j - 1] < sigma[i - 1])
    return math.factorial(N) // prod


def left_covers_down(w):
    """Covers below w in left weak order (Inv-subset order): w -> s_i w."""
    wi = ~w
    return [~(wi.swap(d, d + 1)) for d in wi.zero_indexed_descents()]


def ideal_lambdas(b):
    """Lambda(v) for every v <= b, by bitset DP over the ideal [e,b]."""
    seen = {b}
    frontier = [b]
    while frontier:
        nxt = []
        for w in frontier:
            for x in left_covers_down(w):
                if x not in seen:
                    seen.add(x)
                    nxt.append(x)
        frontier = nxt
    index = {w: i for i, w in enumerate(seen)}
    down = {}
    for w in sorted(seen, key=lambda w: w.inv):
        bits = 1 << index[w]
        for x in left_covers_down(w):
            if x in down:
                bits |= down[x]
        down[w] = bits
    return {w: bin(v).count("1") for w, v in down.items()}


def sweep(n, cap):
    buckets = defaultdict(list)
    for v in Permutation.all_permutations(n):
        buckets[(v.mul_sortable(), v.mul_dominant())].append(v)

    best = None
    best_r = None
    skipped = 0
    max_skipped_lb = 0
    for (a, b), members in buckets.items():
        lb = hook_lambda(b)
        if lb > cap:
            skipped += 1
            max_skipped_lb = max(max_skipped_lb, lb)
            continue
        lam = ideal_lambdas(b)
        mean_lam = sum(lam[v] for v in members) / len(members)
        expo = math.log(lb) / math.log(mean_lam) if mean_lam > 1 else 1.0
        r = lb / mean_lam
        rec = (a, b, len(members), lam[a], lb, mean_lam, r, expo)
        if best is None or expo > best[7]:
            best = rec
        if best_r is None or r > best_r[6]:
            best_r = rec
    return len(buckets), best, best_r, skipped, max_skipped_lb


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 30000

    for n in range(3, 7):
        for v in Permutation.all_permutations(n):
            b = v.mul_dominant()
            if hook_lambda(b) != len(ideal_lambdas(b)):
                raise SystemExit(f"hook formula mismatch at {b}: {hook_lambda(b)} vs {len(ideal_lambdas(b))}")
    print("hook length formula for Lambda(b) verified for n <= 6\n")

    for n in range(5, hi + 1):
        ncls, best, best_r, skipped, mx = sweep(n, cap)
        print(f"n={n}  classes={ncls}  skipped(Lambda(b)>{cap})={skipped} (max skipped Lambda(b)={mx})")
        print(f"   max exponent log Lb/log meanL = {best[7]:.4f}   at {best[0]} .. {best[1]}")
        print(f"      |C|={best[2]} Lambda(a)={best[3]} Lambda(b)={best[4]} mean={best[5]:.2f}")
        print(f"   max ratio Lb/meanL            = {best_r[6]:.4f}   at {best_r[0]} .. {best_r[1]}")
        print(f"      |C|={best_r[2]} Lambda(a)={best_r[3]} Lambda(b)={best_r[4]} mean={best_r[5]:.2f}")
        print(f"      code(b)={list(best_r[1].code)}  code(a)={list(best_r[0].code)}")
        print()

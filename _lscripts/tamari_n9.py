#!/usr/bin/env python3
"""Fast n=9 (and beyond) pass over classes with small Lambda(b).

Works with plain length-n tuples and left weak order covers (swap the values
i, i+1 when i sits to the right of i+1), so no Permutation objects are built in
the inner loops.  Lambda(b) is screened by the hook length formula first, and
the ideal DP is run only for classes under the cap.
"""

import math
import sys
from collections import defaultdict

from schubmult.combinatorics.permutation import Permutation


def hook_lambda(bt):
    n = len(bt)
    inv = [0] * (n + 1)
    for i, x in enumerate(bt):
        inv[x] = i + 1
    sigma = [inv[i] for i in range(1, n + 1)]
    prod = 1
    for i in range(n):
        prod *= 1 + sum(1 for j in range(i) if sigma[j] < sigma[i])
    return math.factorial(n) // prod


def covers_down(wt):
    n = len(wt)
    pos = [0] * (n + 2)
    for i, x in enumerate(wt):
        pos[x] = i
    out = []
    for i in range(1, n):
        pi, pj = pos[i], pos[i + 1]
        if pi > pj:
            lst = list(wt)
            lst[pi], lst[pj] = i + 1, i
            out.append(tuple(lst))
    return out


def ideal_lambdas(bt):
    seen = {bt}
    frontier = [bt]
    while frontier:
        nxt = []
        for w in frontier:
            for x in covers_down(w):
                if x not in seen:
                    seen.add(x)
                    nxt.append(x)
        frontier = nxt
    order = sorted(seen, key=lambda w: sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j]))
    index = {w: i for i, w in enumerate(order)}
    down = {}
    for w in order:
        bits = 1 << index[w]
        for x in covers_down(w):
            d = down.get(x)
            if d is not None:
                bits |= d
        down[w] = bits
    return {w: v.bit_count() for w, v in down.items()}


def sweep(n, cap):
    buckets = defaultdict(list)
    for v in Permutation.all_permutations(n):
        vt = tuple(v[i] for i in range(n))
        a = v.mul_sortable()
        b = v.mul_dominant()
        at = tuple(a[i] for i in range(n))
        bt = tuple(b[i] for i in range(n))
        buckets[(at, bt)].append(vt)

    best_exp = None
    best_r = None
    done = 0
    skipped = 0
    for (at, bt), members in buckets.items():
        lb = hook_lambda(bt)
        if lb > cap:
            skipped += 1
            continue
        done += 1
        lam = ideal_lambdas(bt)
        mean_lam = sum(lam[v] for v in members) / len(members)
        expo = math.log(lb) / math.log(mean_lam) if mean_lam > 1 else 1.0
        r = lb / mean_lam
        rec = (at, bt, len(members), lam[at], lb, mean_lam, r, expo)
        if best_exp is None or expo > best_exp[7]:
            best_exp = rec
        if best_r is None or r > best_r[6]:
            best_r = rec
    return len(buckets), done, skipped, best_exp, best_r


if __name__ == "__main__":
    n = int(sys.argv[1])
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    ncls, done, skipped, be, br = sweep(n, cap)
    print(f"n={n}  classes={ncls}  computed={done}  skipped(Lambda(b)>{cap})={skipped}", flush=True)
    print(f"  max exponent log Lb/log meanL = {be[7]:.4f}")
    print(f"     class {be[0]} .. {be[1]}  |C|={be[2]} La={be[3]} Lb={be[4]} mean={be[5]:.2f}")
    print(f"  max ratio Lb/meanL            = {br[6]:.4f}")
    print(f"     class {br[0]} .. {br[1]}  |C|={br[2]} La={br[3]} Lb={br[4]} mean={br[5]:.2f}")

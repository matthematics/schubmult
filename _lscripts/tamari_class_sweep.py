#!/usr/bin/env python3
"""u-free sweep of the Tamari-class conjectures.

For each v in S_n set a = ddom(v) = v.mul_sortable(), b = mdom(v) = v.mul_dominant(),
group permutations into classes, and test:

  (star)  |C| <= Lambda(a)                     (class no larger than ideal under its min)
  (shadow) Lambda(b) <= poly(n, mean_{v in C} Lambda(v))
  (power)  log|C| / log Lambda(b) < 1 - eps    (power saving in the containment C subset [e,b])

where Lambda(w) = |[e, w]| in the weak order used by Permutation.weak_order_leq.
"""

import sys
from collections import defaultdict

from schubmult.combinatorics.permutation import Permutation


def downset_sizes(n):
    """Lambda(w) for every w in S_n, via bitset dynamic programming up the weak order."""
    perms = Permutation.all_permutations(n)
    index = {w: i for i, w in enumerate(perms)}
    order = sorted(perms, key=lambda w: w.inv)
    down = {}
    for w in order:
        bits = 1 << index[w]
        for d in w.zero_indexed_descents():
            bits |= down[w.swap(d, d + 1)]
        down[w] = bits
    return {w: bin(down[w]).count("1") for w in perms}, index


def classes(n):
    perms = Permutation.all_permutations(n)
    buckets = defaultdict(list)
    for v in perms:
        buckets[(v.mul_sortable(), v.mul_dominant())].append(v)
    return buckets


def check_structure(n, buckets, lam):
    """Sanity checks: a <= v <= b, classes partition S_n, class == interval [a,b]."""
    problems = []
    total = 0
    for (a, b), members in buckets.items():
        total += len(members)
        if a not in members or b not in members:
            problems.append(f"endpoints not in class {a} {b}")
        for v in members:
            if not a.weak_order_leq(v):
                problems.append(f"a !<= v: {a} {v}")
            if not v.weak_order_leq(b):
                problems.append(f"v !<= b: {v} {b}")
    if total != math_factorial(n):
        problems.append(f"class sizes sum to {total}, not {math_factorial(n)}")
    return problems


def math_factorial(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def sweep(n, verbose=False):
    lam, _ = downset_sizes(n)
    buckets = classes(n)
    problems = check_structure(n, buckets, lam)

    rows = []
    for (a, b), members in buckets.items():
        size = len(members)
        la, lb = lam[a], lam[b]
        mean_lam = sum(lam[v] for v in members) / size
        rows.append(
            {
                "a": a,
                "b": b,
                "la_len": a.inv,
                "lb_len": b.inv,
                "diam": b.inv - a.inv,
                "size": size,
                "lam_a": la,
                "lam_b": lb,
                "mean_lam": mean_lam,
                "star": size / la,
                "shadow": lb / mean_lam,
                "power": (0.0 if lb <= 1 else __import__("math").log(size) / __import__("math").log(lb)),
            },
        )

    n_classes = len(rows)
    max_star = max(rows, key=lambda r: r["star"])
    max_shadow = max(rows, key=lambda r: r["shadow"])
    max_power = max(rows, key=lambda r: r["power"])
    max_size = max(rows, key=lambda r: r["size"])

    print(f"n = {n}:  {n_classes} classes, {math_factorial(n)} permutations")
    if problems:
        print("  STRUCTURE PROBLEMS:")
        for p in problems[:10]:
            print(f"    {p}")
    else:
        print("  structure ok (a <= v <= b, classes partition S_n)")
    print(f"  max |C|                = {max_size['size']}  (class {max_size['a']} .. {max_size['b']})")
    print(f"  max |C|/Lambda(a)      = {max_star['star']:.4f}  (star holds iff <= 1)")
    print(f"       at class {max_star['a']} .. {max_star['b']}, |C|={max_star['size']}, Lambda(a)={max_star['lam_a']}")
    print(f"  max Lambda(b)/mean Lam = {max_shadow['shadow']:.4f}")
    print(f"       at class {max_shadow['a']} .. {max_shadow['b']}, Lambda(b)={max_shadow['lam_b']}, mean={max_shadow['mean_lam']:.2f}")
    print(f"  max log|C|/log Lam(b)  = {max_power['power']:.4f}  (power saving iff bounded < 1)")
    print(f"       at class {max_power['a']} .. {max_power['b']}, |C|={max_power['size']}, Lambda(b)={max_power['lam_b']}")
    n_star_fail = sum(1 for r in rows if r["star"] > 1)
    print(f"  classes violating (star): {n_star_fail}")
    print()
    return rows


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    all_rows = {}
    for n in range(2, hi + 1):
        all_rows[n] = sweep(n)

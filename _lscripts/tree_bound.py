#!/usr/bin/env python3
"""Reduce the bounding problem to a statement about binary trees.

For each sylvester (Tamari) class C of S_n, with tree T:

    |C|        = e(T)             (linear extensions of T)  -- product formula
    Lambda(b)  = |[e, mdom v]_L|                            -- product formula
    Lambda(a)  = |[e, ddom v]_L|                            -- ???
    Lbar       = (1/|C|) sum_{v in C} Lambda(v)

Trivially  Lbar >= max( Lambda(a), Lambda(b)/|C| ).

TARGET (sufficient for output-polynomiality):
    max( Lambda(a), Lambda(b)/|C| )  >=  Lambda(b)^(1/c)   for an absolute c.

This script dumps the data so we can (i) find a product formula for Lambda(a)
and (ii) test the target inequality.
"""

import math
import sys
from collections import defaultdict

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")

from count_prods import construct_tree_from_perm, dominant_interval_size  # noqa: E402

from schubmult.combinatorics.permutation import Permutation  # noqa: E402


def shape(t):
    if t is None:
        return None
    return (shape(t.left), shape(t.right))


def size(s):
    if s is None:
        return 0
    return 1 + size(s[0]) + size(s[1])


def ext_count(s):
    """n! / prod of subtree sizes."""
    if s is None:
        return 1
    return math.comb(size(s) - 1, size(s[0])) * ext_count(s[0]) * ext_count(s[1])


def render(s):
    if s is None:
        return "."
    return f"({render(s[0])}{render(s[1])})"


def covers_down(w):
    wi = ~w
    return [~(wi.swap(d, d + 1)) for d in wi.zero_indexed_descents()]


def lam(w):
    """|[e,w]_L| by BFS."""
    seen = {w}
    frontier = [w]
    while frontier:
        nxt = []
        for x in frontier:
            for y in covers_down(x):
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return len(seen)


def main(n, show=12):
    classes = defaultdict(list)
    for v in Permutation.all_permutations(n):
        classes[shape(construct_tree_from_perm(v, n))].append(v)

    lam_cache = {}

    def L(w):
        if w not in lam_cache:
            lam_cache[w] = lam(w)
        return lam_cache[w]

    rows = []
    for s, members in classes.items():
        a = members[0].mul_sortable()
        b = members[0].mul_dominant()
        c = len(members)
        La = L(a)
        Lb = dominant_interval_size(b)
        Lbar = sum(L(x) for x in members) / c
        rows.append((s, c, La, Lb, Lbar))

    # sanity: ext_count matches
    assert all(ext_count(s) == c for s, c, *_ in rows), "ext_count mismatch"

    # target inequality: exponent needed
    def expo(r):
        _, c, La, Lb, _ = r
        best = max(La, Lb / c)
        if Lb <= 1 or best <= 1:
            return 1.0
        return math.log(Lb) / math.log(best)

    def expo_true(r):
        _, c, La, Lb, Lbar = r
        if Lb <= 1 or Lbar <= 1:
            return 1.0
        return math.log(Lb) / math.log(Lbar)

    worst = max(rows, key=expo)
    worst_t = max(rows, key=expo_true)
    print(f"n={n}  classes={len(rows)}")
    print(f"   max exponent vs Lbar          c = {expo_true(worst_t):.4f}   at {render(worst_t[0])}")
    print(f"   max exponent vs max(La,Lb/|C|) c = {expo(worst):.4f}   at {render(worst[0])}")
    print(f"        |C|={worst[1]} La={worst[2]} Lb={worst[3]} Lbar={worst[4]:.2f}")

    print("   shape                     |C|      La       Lb     Lb/|C|    Lbar")
    for s, c, La, Lb, Lbar in sorted(rows, key=expo, reverse=True)[:show]:
        print(f"   {render(s):<24} {c:>6} {La:>7} {Lb:>8} {Lb / c:>9.2f} {Lbar:>9.2f}")
    return rows


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(3, hi + 1):
        main(n)
        print()

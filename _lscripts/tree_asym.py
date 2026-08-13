#!/usr/bin/env python3
"""Asymptotics of |C| = extension_count(T) and the reduced bounding exponent.

(A)  |C| = e(T) = n! / prod_{x in T} h(x),  h(x) = subtree size.
     max_C |C| = n! / F(n),  F(n) = min_T prod h = n * min_{i+j=n-1} F(i)F(j).
     We compute F(n) exactly by DP and extract gamma = lim F(n)^{1/n}.

(B)  d(n) = max_C log|C| / log Lbar,  the exponent that decides output-polynomiality.
     Lambda(v) for all v in S_n by a two-layer bitset DP (low memory).
"""

import math
import sys
from collections import defaultdict
from functools import lru_cache

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")

from count_prods import construct_tree_from_perm, dominant_interval_size  # noqa: E402

from schubmult.combinatorics.permutation import Permutation  # noqa: E402

# ---------------------------------------------------------------- (A)


@lru_cache(maxsize=None)
def F(n):
    """min over binary trees on n nodes of prod of subtree sizes."""
    if n == 0:
        return 1
    return n * min(F(i) * F(n - 1 - i) for i in range(n))


@lru_cache(maxsize=None)
def Fmax(n):
    if n == 0:
        return 1
    return n * max(Fmax(i) * Fmax(n - 1 - i) for i in range(n))


def part_a(hi=200):
    print("(A)  max_C |C| = n!/F(n),  F(n) = min_T prod h(x)")
    print("      n     log2 F(n)   (log2 F)/n     n!/F(n) = max|C|")
    for n in [7, 15, 31, 63, 100, 150, hi]:
        lf = math.log2(F(n))
        print(f"   {n:>5}   {lf:>10.3f}   {lf / n:>10.5f}    ~ n!/{2 ** (lf / n):.4f}^n")
    g = 2 ** (math.log2(F(hi)) / hi)
    print(f"\n   gamma = lim F(n)^(1/n) ~ {g:.5f}      (min |C| = 1, at the comb)")
    print(f"   so  max_C |C| = n! / {g:.4f}^(n(1+o(1)))   and  sum_C |C| = n!, #C = Catalan ~ 4^n")


# ---------------------------------------------------------------- (B)


def shape(t):
    if t is None:
        return None
    return (shape(t.left), shape(t.right))


def render(s):
    return "." if s is None else f"({render(s[0])}{render(s[1])})"


def all_lambdas(n):
    """Lambda(v)=|[e,v]_L| for all v in S_n, keeping only two length layers of bitsets."""
    perms = list(Permutation.all_permutations(n))
    idx = {p: i for i, p in enumerate(perms)}
    layers = defaultdict(list)
    for p in perms:
        layers[p.inv].append(p)
    lam = {}
    prev = {}
    for k in sorted(layers):
        cur = {}
        for w in layers[k]:
            wi = ~w
            bits = 1 << idx[w]
            for d in wi.zero_indexed_descents():
                bits |= prev[~(wi.swap(d, d + 1))]
            cur[w] = bits
            lam[w] = bin(bits).count("1")
        prev = cur
    return lam


def part_b(hi):
    print("\n(B)  reduced bounding exponent  d(n) = max_C log|C| / log Lbar")
    print("      n   #classes      d(n)     c(n)=logLb/logLbar   argmax tree")
    prev_d = None
    for n in range(3, hi + 1):
        lam = all_lambdas(n)
        classes = defaultdict(list)
        for v in Permutation.all_permutations(n):
            classes[shape(construct_tree_from_perm(v, n))].append(v)
        best_d, best_s, best_c = 0.0, None, 0.0
        for s, mem in classes.items():
            c = len(mem)
            Lbar = sum(lam[x] for x in mem) / c
            if Lbar > 1:
                if c > 1 and math.log(c) / math.log(Lbar) > best_d:
                    best_d = math.log(c) / math.log(Lbar)
                    best_s = s
                Lb = dominant_interval_size(mem[0].mul_dominant())
                if Lb > 1:
                    best_c = max(best_c, math.log(Lb) / math.log(Lbar))
        inc = "" if prev_d is None else f"  (+{best_d - prev_d:.4f})"
        prev_d = best_d
        print(f"   {n:>4} {len(classes):>9}   {best_d:.4f}{inc:<12} {best_c:.4f}       {render(best_s)}")


if __name__ == "__main__":
    part_a()
    part_b(int(sys.argv[1]) if len(sys.argv) > 1 else 7)

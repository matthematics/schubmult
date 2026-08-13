#!/usr/bin/env python3
"""THE bounding problem, reduced to binary trees.

Chain (paper):
    W          = O(n^4 Lambda(b)^2 N(u,b))      and Lambda(b) <= N(u,b)
               = O(n^4 (N^up)^3)
    N^up       = N(u,b) <= |C| * Ntilde         (b in C, Ntilde is the average)
    Ntilde     >= Lbar   (termwise N(u,v') >= Lambda(v'))
  =>  W = O(n^4 |C|^3 Ntilde^3)

So output-polynomiality  <==>   |C| <= poly(n) * Lbar^d   for absolute d.

Here |C| = e(T) is an exact O(n) product over the tree.  This script measures
    d(C) = log|C| / log Lbar
and also tests whether S(T) = sum_{v in C} Lambda(v) obeys a tree recursion.
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
    return 0 if s is None else 1 + size(s[0]) + size(s[1])


def ext_count(s):
    if s is None:
        return 1
    return math.comb(size(s) - 1, size(s[0])) * ext_count(s[0]) * ext_count(s[1])


def render(s):
    return "." if s is None else f"({render(s[0])}{render(s[1])})"


def all_lambdas(n):
    """Lambda(v) = |[e,v]_L| for every v in S_n, by ideal-union DP over length layers."""
    perms = list(Permutation.all_permutations(n))
    idx = {p: i for i, p in enumerate(perms)}
    perms.sort(key=lambda p: p.inv)
    ideals = {}
    lam = {}
    for w in perms:
        wi = ~w
        bits = 1 << idx[w]
        for d in wi.zero_indexed_descents():
            bits |= ideals[~(wi.swap(d, d + 1))]
        ideals[w] = bits
        lam[w] = bin(bits).count("1")
    return lam


def main(n):
    lam = all_lambdas(n)
    classes = defaultdict(list)
    for v in Permutation.all_permutations(n):
        classes[shape(construct_tree_from_perm(v, n))].append(v)

    rows = []
    for s, members in classes.items():
        c = len(members)
        assert ext_count(s) == c
        S = sum(lam[x] for x in members)
        Lbar = S / c
        Lb = dominant_interval_size(members[0].mul_dominant())
        La = lam[members[0].mul_sortable()]
        rows.append((s, c, S, Lbar, La, Lb))

    def d_of(r):
        _, c, _, Lbar, _, _ = r
        if c <= 1:
            return 0.0
        return math.log(c) / math.log(Lbar) if Lbar > 1 else float("inf")

    def c_of(r):
        _, _, _, Lbar, _, Lb = r
        return math.log(Lb) / math.log(Lbar) if Lbar > 1 and Lb > 1 else 1.0

    wd = max(rows, key=d_of)
    wc = max(rows, key=c_of)
    print(f"n={n}  classes={len(rows)}")
    print(f"   max d = log|C|/log Lbar      = {d_of(wd):.4f}   |C|={wd[1]} Lbar={wd[3]:.2f}  {render(wd[0])}")
    print(f"   max c = log Lb /log Lbar     = {c_of(wc):.4f}   Lb={wc[5]} Lbar={wc[3]:.2f}  {render(wc[0])}")
    print(f"   max |C|/Lbar                 = {max(r[1] / r[3] for r in rows):.3f}")
    print(f"   max (|C|*Lb)/Lbar^2          = {max(r[1] * r[5] / r[3] ** 2 for r in rows):.3f}")
    return {s: (c, S) for s, c, S, *_ in rows}


def test_recursion(data):
    """Is S(T) = sum_{v in C_T} Lambda(v) determined by S(left), S(right)?

    Guess: S(T) = binom(n-1, |L|) * S(L) * S(R) * g(n, |L|) for some simple g.
    Print the residual ratio.
    """
    print("\n  testing tree recursion for S(T) = sum Lambda(v):")
    shown = 0
    for s, (c, S) in sorted(data.items(), key=lambda kv: -kv[1][1]):
        if s is None or shown >= 10:
            continue
        left, right = s
        dl = data.get(left)
        dr = data.get(right)
        if dl is None or dr is None:
            continue
        n = size(s)
        base = math.comb(n - 1, size(left)) * dl[1] * dr[1]
        print(f"      {render(s):<22} S={S:<7} base={base:<9} S/base={S / base:.5f}")
        shown += 1


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    last = None
    for n in range(3, hi + 1):
        last = main(n)
        print()
    # recursion test needs all n simultaneously; rebuild cumulative table
    cum = {}
    for n in range(1, hi + 1):
        lam = all_lambdas(n)
        cl = defaultdict(list)
        for v in Permutation.all_permutations(n):
            cl[shape(construct_tree_from_perm(v, n))].append(v)
        for s, mem in cl.items():
            cum[s] = (len(mem), sum(lam[x] for x in mem))
    test_recursion(cum)

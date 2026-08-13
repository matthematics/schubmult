#!/usr/bin/env python3
"""Validate the tree formulas and hunt for Lambda(b) in terms of the tree.

Claims under test, for v in S_n with T = construct_tree_from_perm(v, n):

  (1) T determines the Tamari (sylvester) class of v.
  (2) T.extension_count() == |[ddom v, mdom v]_L|.
  (3) dominant_interval_size(mdom v) == Lambda(mdom v) = |[e, mdom v]_L|.

and then we look for a product formula expressing Lambda(b) over the tree.
"""

import sys
from collections import defaultdict

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")

from count_prods import TreeNode, construct_tree_from_perm, dominant_interval_size  # noqa: E402

from schubmult.combinatorics.permutation import Permutation  # noqa: E402


def shape(t):
    if t is None:
        return None
    return (shape(t.left), shape(t.right))


def covers_down(w):
    wi = ~w
    return [~(wi.swap(d, d + 1)) for d in wi.zero_indexed_descents()]


def ideal(b):
    seen = {b}
    frontier = [b]
    while frontier:
        nxt = []
        for w in frontier:
            for x in covers_down(w):
                if x not in seen:
                    seen.add(x)
                    nxt.append(x)
        frontier = nxt
    return seen


def hooks(t):
    """Multiset of subtree sizes of t."""
    if t is None:
        return []
    return [t.total_size(), *hooks(t.left), *hooks(t.right)]


def main(n):
    perms = Permutation.all_permutations(n)

    by_class = defaultdict(list)
    by_shape = defaultdict(list)
    for v in perms:
        by_class[(v.mul_sortable(), v.mul_dominant())].append(v)
        by_shape[shape(construct_tree_from_perm(v, n))].append(v)

    same = set(map(lambda L: tuple(sorted(L, key=str)), by_class.values())) == set(map(lambda L: tuple(sorted(L, key=str)), by_shape.values()))
    print(f"n={n}: {len(by_class)} classes, {len(by_shape)} tree shapes; partitions agree: {same}")

    bad_ext = 0
    bad_lam = 0
    rows = []
    for (a, b), members in by_class.items():
        t = construct_tree_from_perm(members[0], n)
        ec = t.extension_count()
        if ec != len(members):
            bad_ext += 1
            if bad_ext <= 3:
                print(f"   extension_count {ec} != |C| {len(members)}  at {a} .. {b}")
        lam_b_true = len(ideal(b))
        lam_b_form = dominant_interval_size(b)
        if lam_b_true != lam_b_form:
            bad_lam += 1
            if bad_lam <= 3:
                print(f"   dominant_interval_size {lam_b_form} != Lambda(b) {lam_b_true}  at b={b}")
        rows.append((a, b, len(members), ec, lam_b_true, sorted(hooks(t), reverse=True)))

    print(f"   extension_count mismatches: {bad_ext}   Lambda(b) mismatches: {bad_lam}")

    print("   sample (|C|, Lambda(b), hooks of T):")
    for r in sorted(rows, key=lambda r: -r[4])[:8]:
        print(f"      |C|={r[2]:<5} Lambda(b)={r[4]:<6} hooks={r[5]}  b={r[1]}")
    return rows


if __name__ == "__main__":
    for n in range(3, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 7):
        main(n)
        print()

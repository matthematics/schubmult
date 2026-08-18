"""Test the mechanism behind the v1 (independent) right-ideal argument.

v1 proves Gamma is a right ideal by working one commutation move at a time and
pairing terms in the fiber by their word suffixes.  The natural inductive form
of that claim is:

  (LIFT-MOVE)  Let R_1, R_2 have height m, the same permutation u, and let
  word(R_2) be obtained from word(R_1) by a SINGLE adjacent commutation that
  preserves the forest class.  For w with w \downvar{m+1} u (equal length) let
  B_w, B'_w be the unique elements of Z^{-1}(R_1), Z^{-1}(R_2) in RC_{m+1}(w)^0.
  Then word(B'_w) is obtained from word(B_w) by a single adjacent commutation
  (which then automatically preserves the forest class).

We test:
  (a) whether the fibers really are indexed by the same permutation set,
  (b) LIFT-MOVE,
  (c) the weaker statement actually needed: B_w =_forest B'_w,
  (d) the "commutation distance" between word(B_w) and word(B'_w).
"""

import sys
from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.combinatorics.rc_graph import RCGraph


def pair_label(word):
    counts = {}
    out = []
    for a in word:
        a = int(a)
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def finv(word):
    return omega_insertion(pair_label(tuple(reversed(word))))[0]


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def comm_neighbors(word):
    out = []
    for i in range(len(word) - 1):
        if abs(word[i] - word[i + 1]) >= 2:
            out.append(word[:i] + (word[i + 1], word[i]) + word[i + 2 :])
    return out


def comm_distance(x, y, cap=4):
    """number of adjacent commutations needed to turn word x into word y"""
    if x == y:
        return 0
    frontier = {x}
    seen = {x}
    for d in range(1, cap + 1):
        nxt = set()
        for u in frontier:
            for v in comm_neighbors(u):
                if v == y:
                    return d
                if v not in seen:
                    seen.add(v)
                    nxt.add(v)
        frontier = nxt
        if not frontier:
            return None
    return None


def main(N=5):
    idx_tot = idx_bad = 0
    lift_tot = lift_bad = 0
    need_tot = need_bad = 0
    disthist = defaultdict(int)
    examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            u = Permutation(list(pw))
            d = maxd(u)
            if d == 0:
                continue
            for m in range(max(d, 1), n + 1):
                if (u, m) in seen:
                    continue
                seen.add((u, m))
                graphs = rc0(u, m + 1)
                # fiber map, organized by the downstairs graph
                fiber = defaultdict(dict)   # R -> {perm w : B}
                for B in graphs:
                    fiber[B.zero_out_last_row()][B.perm] = B
                # group downstairs graphs of height m by forest class
                byclass = defaultdict(list)
                for R in fiber:
                    byclass[R.forest_invariant].append(R)
                for phi, Rs in byclass.items():
                    for i in range(len(Rs)):
                        for j in range(len(Rs)):
                            if i == j:
                                continue
                            R1, R2 = Rs[i], Rs[j]
                            # only single commutation moves apart
                            if tuple(R2.perm_word) not in comm_neighbors(tuple(R1.perm_word)):
                                continue
                            k1, k2 = fiber[R1], fiber[R2]
                            idx_tot += 1
                            if set(k1) != set(k2):
                                idx_bad += 1
                                continue
                            for w in k1:
                                B, Bp = k1[w], k2[w]
                                x, y = tuple(B.perm_word), tuple(Bp.perm_word)
                                need_tot += 1
                                if B.forest_invariant != Bp.forest_invariant:
                                    need_bad += 1
                                lift_tot += 1
                                dist = comm_distance(x, y)
                                disthist[dist] += 1
                                if dist != 1:
                                    lift_bad += 1
                                    if len(examples) < 5:
                                        examples.append((tuple(u[:n]), m, tuple(w[: n + 1]), x, y, dist))
    print(f"fibers indexed by same permutation set: {idx_tot - idx_bad}/{idx_tot} (bad {idx_bad})")
    print(f"LIFT-MOVE (one commutation lifts to one commutation): {lift_tot - lift_bad}/{lift_tot} (bad {lift_bad})")
    print(f"   commutation-distance histogram (None = >4 or not commutation-equivalent): {dict(sorted(disthist.items(), key=lambda t: (t[0] is None, t[0])))}")
    print(f"needed statement B =_forest B': {need_tot - need_bad}/{need_tot} (bad {need_bad})")
    for e in examples:
        print("   LIFT-MOVE fails:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

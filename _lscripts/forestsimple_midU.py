"""Isolate the problem case a = b+2, a' = a, b' = b+1 of eq:forestsimple_positions.

Among instances with downstairs gap |a-b| = 2, compare the forest-equivalent ones
(where the claim holds) with the failures, looking for a distinguishing feature of
the prefix U: in particular whether the intermediate value v+1 occurs in U.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def swap_positions(w1, w2):
    if len(w1) != len(w2) or w1 == w2:
        return None
    diff = [i for i in range(len(w1)) if w1[i] != w2[i]]
    if len(diff) != 2:
        return None
    p, q = diff
    if q != p + 1:
        return None
    if w1[p] == w2[q] and w1[q] == w2[p]:
        return p
    return None


def main(n):
    stats = Counter()

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        maxd = max(w.descents()) + 1
        for h in range(maxd, n + 1):
            if maxd <= h - 1:
                continue
            graphs = [B for B in RCGraph.all_rc_graphs(w, h) if len(B[-1]) == 0]
            images = {}
            for B in graphs:
                try:
                    images[B] = B.zero_out_last_row()
                except Exception:
                    pass
            for B1, B2 in itertools.permutations(images, 2):
                R1, R2 = images[B1], images[B2]
                if R1.perm != R2.perm:
                    continue
                W = R1.perm_word
                p = swap_positions(W, R2.perm_word)
                if p is None:
                    continue
                a, b = W[p], W[p + 1]
                if abs(a - b) != 2:
                    continue
                v = min(a, b)  # values are v and v+2
                U = W[:p]
                feq = R1.forest_invariant == R2.forest_invariant
                gap_up = abs(B1.perm_word[p] - B1.perm_word[p + 1])
                stats[(
                    "feq" if feq else "noteq",
                    f"gap_up={gap_up}",
                    "mid_in_U" if (v + 1) in U else "mid_not_in_U",
                )] += 1

    print(f"n = {n}   (instances with downstairs gap |a-b| = 2; mid = v+1)")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

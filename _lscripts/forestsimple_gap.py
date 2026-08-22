"""Probe the mechanism behind eq:forestsimple_positions.

For each nontrivial instance, record the integer gap |a-b| of the swapped letters
downstairs and upstairs, split by forest equivalence and by success. The proposed
mechanism: two adjacent letters with |a-b| >= 2 give disjoint crossings, so
swapping them is an isotopy of the wiring diagram and the Little bump acts
identically; the claim can only fail if the bump brings the letters to distance 1.
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
    gap_stats = Counter()
    fail_gap = Counter()

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
                p = swap_positions(R1.perm_word, R2.perm_word)
                if p is None:
                    continue
                feq = R1.forest_invariant == R2.forest_invariant
                ok = swap_positions(B1.perm_word, B2.perm_word) == p
                gap_down = abs(R1.perm_word[p] - R1.perm_word[p + 1])
                gap_up = abs(B1.perm_word[p] - B1.perm_word[p + 1])
                tag = "feq" if feq else "noteq"
                gap_stats[(tag, "ok" if ok else "FAIL", f"gap_up={gap_up}")] += 1
                if not ok:
                    fail_gap[(tag, f"gap_down={gap_down}")] += 1

    print(f"n = {n}  (nontrivial bump cases; gap_up = |a-b| in the upstairs word)")
    for k in sorted(gap_stats, key=str):
        print(f"  {k}: {gap_stats[k]}")
    print("\ndownstairs gap among failures:")
    for k in sorted(fail_gap, key=str):
        print(f"  {k}: {fail_gap[k]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

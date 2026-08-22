"""Classify the eq:forestsimple_positions instances by swap type.

For each pair B1,B2 in RC_{m+1}(w)^0 whose images under Z differ by an adjacent
transposition at positions p,p+1, record |a-b| downstairs (commutation vs braid),
whether the images are forest equivalent, and whether the positions are preserved
upstairs.
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
    braid_feq = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        maxd = max(w.descents()) + 1
        for h in range(maxd, n + 1):
            if maxd <= h - 1:
                continue  # only the genuine bump case
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
                a, b = R1.perm_word[p], R1.perm_word[p + 1]
                kind = "comm" if abs(a - b) >= 2 else "braid"
                feq = R1.forest_invariant == R2.forest_invariant
                ok = swap_positions(B1.perm_word, B2.perm_word) == p
                stats[(kind, "feq" if feq else "noteq", "ok" if ok else "FAIL")] += 1
                if kind == "braid" and feq and len(braid_feq) < 3:
                    braid_feq.append((w, h, R1, R2))

    print(f"n = {n}  (nontrivial bump cases only)")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if braid_feq:
        print("\nbraid + forest-equivalent examples:")
        for w, h, R1, R2 in braid_feq:
            print(f"  w={w} h={h} {R1.perm_word} vs {R2.perm_word}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

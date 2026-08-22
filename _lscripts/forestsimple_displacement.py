"""Check the per-letter displacement of Z on the two swapped positions.

Little's Lemma 5 gives uniqueness of decremented indices within a single bump,
but Z iterates bumps, so a priori a letter could move several times. This records
the actual displacements a' - a and b' - b at the swapped positions.
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
    allpos = Counter()

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        maxd = max(w.descents()) + 1
        for h in range(maxd, n + 1):
            if maxd <= h - 1:
                continue
            for B in RCGraph.all_rc_graphs(w, h):
                if len(B[-1]) != 0:
                    continue
                try:
                    R = B.zero_out_last_row()
                except Exception:
                    continue
                # displacement over every position, not just swapped ones
                for i, (u, d) in enumerate(zip(B.perm_word, R.perm_word)):
                    allpos[u - d] += 1

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
                da = B1.perm_word[p] - R1.perm_word[p]
                db = B1.perm_word[p + 1] - R1.perm_word[p + 1]
                feq = R1.forest_invariant == R2.forest_invariant
                stats[("feq" if feq else "noteq", f"(a'-a,b'-b)=({da},{db})")] += 1

    print(f"n = {n}")
    print("displacement u - d over ALL positions of ALL graphs:")
    for k in sorted(allpos):
        print(f"  {k}: {allpos[k]}")
    print("\ndisplacements at the swapped positions:")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

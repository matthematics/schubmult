"""How much of the prefix does the bump actually touch?

For the forest-equivalent instances, the upstairs prefix U' and the downstairs
prefix U (the letters before position p) differ only at the positions the bump
decrements. If the bump never touches the prefix, then U' = U verbatim and
rl(U') = rl(U), which would be a usable base case for an induction on the number
of decremented positions lying in the prefix.
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


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            if md <= h - 1:
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
                if p is None or R1.forest_invariant != R2.forest_invariant:
                    continue
                up, down = B1.perm_word, R1.perm_word
                # positions strictly before p that the bump decremented
                touched = [i for i in range(p) if up[i] != down[i]]
                gap = abs(down[p] - down[p + 1])
                stats[f"prefix decrements = {len(touched)}"] += 1
                stats[f"  gap={gap}, prefix decrements = {len(touched)}"] += 1
                # also: how many letters of the whole word were decremented
                total = sum(1 for i in range(len(up)) if up[i] != down[i])
                stats[f"    total decrements = {total}"] += 1

    print(f"n = {n}  (forest-equivalent nontrivial instances)")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

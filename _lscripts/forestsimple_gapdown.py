"""Does forest equivalence force the downstairs gap |a-b| >= 3?

If so, then since a Little bump alters each letter by at most one, the upstairs
gap satisfies |a'-b'| >= |a-b| - 1 >= 2, which is eq:forestsimple_positions.
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
    examples = []

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
                gap_down = abs(R1.perm_word[p] - R1.perm_word[p + 1])
                gap_up = abs(B1.perm_word[p] - B1.perm_word[p + 1])
                tag = "feq" if feq else "noteq"
                stats[(tag, f"gap_down={gap_down}", f"gap_up={gap_up}")] += 1
                if feq and gap_down == 2 and len(examples) < 5:
                    examples.append((w, h, B1.perm_word, R1.perm_word, p))

    print(f"n = {n}")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if examples:
        print("\nforest-equivalent instances with gap_down = 2:")
        for w, h, wb, wr, p in examples:
            print(f"  w={w} h={h} p={p} word(B1)={wb} word(R1)={wr}")
    else:
        print("\nNO forest-equivalent instance has gap_down = 2.")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

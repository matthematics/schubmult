"""Probe stronger forms of eq:forestsimple_positions on the forest-equivalent instances.

Records, for each forest-equivalent nontrivial instance: whether the upstairs words
agree away from positions p, p+1; whether the letters at those positions are the
same pair upstairs as downstairs; and whether the bumps start at the same position.
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
                p = swap_positions(R1.perm_word, R2.perm_word)
                if p is None or R1.forest_invariant != R2.forest_invariant:
                    continue
                stats["total_feq"] += 1

                u1, u2 = B1.perm_word, B2.perm_word
                # agreement away from p, p+1
                away = all(u1[i] == u2[i] for i in range(len(u1)) if i not in (p, p + 1))
                stats["agree_away_from_p"] += bool(away)
                # same unordered letter pair at p,p+1 upstairs and downstairs
                same_pair = {u1[p], u1[p + 1]} == {R1.perm_word[p], R1.perm_word[p + 1]}
                stats["letters_unchanged_at_p"] += bool(same_pair)
                # upstairs is also a commutation
                a, b = u1[p], u1[p + 1]
                stats["upstairs_commutation"] += bool(abs(a - b) >= 2)
                # forest equivalent upstairs
                stats["upstairs_feq"] += bool(B1.forest_invariant == B2.forest_invariant)
                # whole word equal except swap (implied by away + swap)
                stats["exact_swap_upstairs"] += bool(swap_positions(u1, u2) == p)

    print(f"n = {n}  (forest-equivalent nontrivial bump instances)")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

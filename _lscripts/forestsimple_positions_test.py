"""Empirically test eq:forestsimple_positions of forest_lr_rule.tex.

Claim: if B1, B2 in RC_{m+1}(w)^0 have Z(B1)=R1, Z(B2)=R2 with word(R2) obtained
from word(R1) by interchanging positions p, p+1, then word(B2) is word(B1) with
the entries in the same positions p, p+1 interchanged.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def swap_positions(w1, w2):
    """Return p if w2 is w1 with entries at positions p, p+1 interchanged, else None."""
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
    failures = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        maxd = max(w.descents()) + 1
        for h in range(maxd, n + 1):
            # B in RC_h(w) with empty last row
            graphs = [B for B in RCGraph.all_rc_graphs(w, h) if len(B[-1]) == 0]
            if len(graphs) < 2:
                continue
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
                p_down = swap_positions(R1.perm_word, R2.perm_word)
                if p_down is None:
                    continue
                forest_eq = R1.forest_invariant == R2.forest_invariant
                trivial = maxd <= h - 1  # Z just deletes the empty last row
                key = ("trivial" if trivial else "bump", "feq" if forest_eq else "noteq")
                stats[key] += 1

                p_up = swap_positions(B1.perm_word, B2.perm_word)
                if p_up == p_down:
                    stats[(*key, "same_positions")] += 1
                else:
                    stats[(*key, "FAIL")] += 1
                    if forest_eq and not trivial and len(failures) < 5:
                        failures.append((w, h, B1, B2, R1, R2, p_down, p_up))

    print(f"n = {n}")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if failures:
        print("\nsample failures (forest-equivalent, nontrivial):")
        for w, h, B1, B2, R1, R2, p_down, p_up in failures:
            print(f"  w={w} h={h} p_down={p_down} p_up={p_up}")
            print(f"    word(B1)={B1.perm_word} word(B2)={B2.perm_word}")
            print(f"    word(R1)={R1.perm_word} word(R2)={R2.perm_word}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

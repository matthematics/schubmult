"""Is the separator set itself monotone, restricted to the interval?

Global containment of rootlists fails, but the only thing eq:sepmono needs is the
elements lying strictly between the two letters. Tests, per case, whether the
downstairs separators inject into the upstairs ones.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph

LO, HI = -4, 40


def rootlist(vals):
    rl = sorted(range(LO, HI))
    for a in vals:
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


def seps(a, b, rl):
    lo, hi = (a, b) if a < b else (b, a)
    return {r for r in rl if lo < r < hi}


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()
    ex = []

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
                d1, d2 = R1.perm_word, R2.perm_word
                diff = [i for i in range(len(d1)) if d1[i] != d2[i]]
                if len(diff) != 2 or diff[1] != diff[0] + 1:
                    continue
                p = diff[0]
                if not (d1[p] == d2[p + 1] and d1[p + 1] == d2[p]):
                    continue
                if R1.forest_invariant != R2.forest_invariant:
                    continue

                up, down = B1.perm_word, d1
                a, b = down[p], down[p + 1]
                ap, bp = up[p], up[p + 1]
                S = seps(a, b, rootlist(list(reversed(down[p + 2 :]))))
                Sp = seps(ap, bp, rootlist(list(reversed(up[p + 2 :]))))
                case = f"da={ap-a},db={bp-b}"
                stats[(case, "separators contained" if S <= Sp else "NOT contained")] += 1
                if not S <= Sp and len(ex) < 6:
                    ex.append((up, down, p, sorted(S), sorted(Sp)))

    print(f"n = {n}  (forest-equivalent instances)")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if ex:
        print("\nnon-containment of separators:")
        for up, down, p, S, Sp in ex:
            print(f"  up={up} down={down} p={p}: down {S} vs up {Sp}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

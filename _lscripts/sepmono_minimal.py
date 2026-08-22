"""Isolate the minimal case of eq:sepmono that actually needs an argument.

For each forest-equivalent instance, record how the two endpoints move under the
bump and whether the segment after the pair is altered at all. Cases where the
segment is unchanged and both endpoints are fixed are trivial; the residue is the
statement that has to be proved.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph

LO, HI = -4, 40


def rootlist(vals):
    """Untagged rootlist: insert bare values in the given order."""
    rl = sorted(range(LO, HI))
    for a in vals:
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


def sep(a, b, rl):
    lo, hi = (a, b) if a < b else (b, a)
    return sum(1 for r in rl if lo < r < hi)


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()
    hard = []

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
                U = list(reversed(down[p + 2 :]))
                Up = list(reversed(up[p + 2 :]))

                da, db = ap - a, bp - b
                tail_same = U == Up
                s_dn = sep(a, b, rootlist(U))
                s_up = sep(ap, bp, rootlist(Up))

                key = (f"da={da},db={db}", "tail same" if tail_same else "tail changed")
                stats[key] += 1
                if s_up < s_dn:
                    stats[("COUNT DROPS", *key)] += 1
                if not (tail_same and da == 0 and db == 0):
                    stats["nontrivial"] += 1
                    if s_dn >= 2 and s_up == 2 and len(hard) < 8:
                        hard.append((up, down, p, da, db, tail_same, s_dn, s_up))
                else:
                    stats["trivial (tail same, endpoints fixed)"] += 1

    print(f"n = {n}  (forest-equivalent instances)")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if hard:
        print("\ntight cases (upstairs count exactly 2):")
        for up, down, p, da, db, ts, sd, su in hard:
            print(f"  up={up} down={down} p={p} da={da} db={db} tail_same={ts} sep {sd}->{su}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

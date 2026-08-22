"""Does the rootlist grow going upward through a bump?

Tests whether rl(U') contains rl(U) when U' bumps to U, which would give the
monotonicity of the separator count immediately. Compares on the window spanned
by the letters actually involved.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph

LO, HI = -4, 40


def injectify(word):
    seen = Counter()
    out = []
    for v in word:
        seen[v] += 1
        out.append((v, seen[v]))
    return out


def rootlist(letters):
    rl = sorted((i, 0) for i in range(LO, HI))
    for a in letters:
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


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
                rl_dn = set(rootlist(list(reversed(injectify(down)[p + 2 :]))))
                rl_up = set(rootlist(list(reversed(injectify(up)[p + 2 :]))))

                stats["rl(U) subset of rl(U')" if rl_dn <= rl_up else "NOT a subset"] += 1
                stats["rl(U') subset of rl(U)" if rl_up <= rl_dn else "reverse not a subset"] += 1
                stats["rootlists equal" if rl_up == rl_dn else "rootlists differ"] += 1
                if not (rl_dn <= rl_up) and len(ex) < 4:
                    lost = sorted(rl_dn - rl_up)
                    gained = sorted(rl_up - rl_dn)
                    ex.append((up, down, p, lost, gained))

    print(f"n = {n}  (forest-equivalent instances)")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if ex:
        print("\nnon-containment examples:")
        for up, down, p, lost, gained in ex:
            print(f"  up={up} down={down} p={p}")
            print(f"    in rl(U) not rl(U'): {lost}")
            print(f"    in rl(U') not rl(U): {gained}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

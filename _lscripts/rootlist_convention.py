"""Determine the correct conventions for the rootlist simulation.

Tries the four combinations of (insertion order) x (tag assignment order) and
reports which reproduces forest equivalence via the Proposition 5.8 criterion.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph

LO, HI = -4, 40


def injectify(word, reverse_tags):
    seen = Counter()
    idx = range(len(word) - 1, -1, -1) if reverse_tags else range(len(word))
    out = [None] * len(word)
    for i in idx:
        v = word[i]
        seen[v] += 1
        out[i] = (v, seen[v])
    return out


def rootlist(letters):
    rl = sorted((i, 0) for i in range(LO, HI))
    for a in letters:
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        if not below or not above:
            raise RuntimeError("window too small")
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


def separated(a, b, rl):
    lo, hi = (a, b) if a < b else (b, a)
    return sum(1 for r in rl if lo < r < hi) >= 2


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    combos = [(ro, rt) for ro in (False, True) for rt in (False, True)]
    stats = {c: Counter() for c in combos}

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            graphs = [B for B in RCGraph.all_rc_graphs(w, h) if len(B[-1]) == 0]
            imgs = []
            for B in graphs:
                try:
                    imgs.append(B.zero_out_last_row())
                except Exception:
                    pass
            for R1, R2 in itertools.permutations(imgs, 2):
                if R1.perm != R2.perm:
                    continue
                d1, d2 = R1.perm_word, R2.perm_word
                diff = [i for i in range(len(d1)) if d1[i] != d2[i]]
                if len(diff) != 2 or diff[1] != diff[0] + 1:
                    continue
                p = diff[0]
                if not (d1[p] == d2[p + 1] and d1[p + 1] == d2[p]):
                    continue
                actual = R1.forest_invariant == R2.forest_invariant
                for ro, rt in combos:
                    inj = injectify(d1, rt)
                    if ro:
                        # NT prefix = the part inserted first = suffix of our word
                        pre = list(reversed(inj[p + 2 :]))
                    else:
                        pre = inj[:p]
                    pred = separated(inj[p], inj[p + 1], rootlist(pre))
                    stats[(ro, rt)]["agree" if pred == actual else "MISMATCH"] += 1

    print(f"n = {n}")
    for c in combos:
        lbl = f"insert={'reverse' if c[0] else 'forward'}, tags={'reverse' if c[1] else 'forward'}"
        print(f"  {lbl}: {dict(stats[c])}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

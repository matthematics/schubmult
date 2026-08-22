"""Test the local window lemma for the effect of a bump on the rootlist.

Claim: if W' bumps to W, decrementing the letters at positions D (values v -> v-1),
then rl(W') and rl(W) agree outside the union of the intervals [v-1, v+1] over the
decremented values v.

Also tested, as the sharper form: outside that union the two rootlists agree
elementwise, and the bare elements freed/consumed are exactly v+1 and v-1.
"""

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


def rootlist(word):
    rl = sorted((i, 0) for i in range(LO, HI))
    for a in reversed(injectify(word)):
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()
    bad = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            if md != h:
                continue
            for B in RCGraph.all_rc_graphs(w, h):
                if len(B[-1]) != 0:
                    continue
                nxt = B.little_bump_desc()
                if nxt == B:
                    continue
                W, Wb = B.perm_word, nxt.perm_word
                if len(W) != len(Wb):
                    continue
                dec = [i for i in range(len(W)) if W[i] != Wb[i]]
                if not dec:
                    stats["no decrement"] += 1
                    continue
                if any(W[i] - Wb[i] != 1 for i in dec):
                    stats["NOT a unit decrement"] += 1
                    continue

                rl_before, rl_after = rootlist(W), rootlist(Wb)
                vals = {W[i] for i in dec}
                windowed = set()
                for v in vals:
                    windowed |= {v - 1, v, v + 1}

                sym = set(rl_before) ^ set(rl_after)
                if all(r[0] in windowed for r in sym):
                    stats["agree outside the windows"] += 1
                else:
                    stats["ESCAPES the windows"] += 1
                    if len(bad) < 6:
                        outside = sorted(r for r in sym if r[0] not in windowed)
                        bad.append((W, Wb, sorted(vals), outside))

                # sharper: for a single decrement, is v+1 freed and v-1 consumed?
                if len(dec) == 1:
                    v = W[dec[0]]
                    freed = (v + 1, 0) in set(rl_after) - set(rl_before)
                    eaten = (v - 1, 0) in set(rl_before) - set(rl_after)
                    stats[f"single: freed(v+1)={freed}, consumed(v-1)={eaten}"] += 1

    print(f"n = {n}  (single bumps)")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if bad:
        print("\nescapes:")
        for W, Wb, vals, outside in bad:
            print(f"  {W} -> {Wb}, decremented values {vals}")
            print(f"    symmetric difference outside windows: {outside}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

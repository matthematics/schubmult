"""What does a single Little bump do to the forest of the whole word?

Uses the validated rootlist simulation (NT insertion order = reverse of ours,
Lemma 5.7) to compare rl(W) with rl(W-bumped), and schubmult's forest data to
compare the supports and the P-symbols.
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
    """rl of the whole word, inserted in NT order (reverse of ours)."""
    rl = sorted((i, 0) for i in range(LO, HI))
    for a in reversed(injectify(word)):
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


def window(rl, lo=0, hi=12):
    return [r for r in rl if lo <= r[0] <= hi]


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()
    examples = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            if md != h:
                continue  # only genuine bumps
            for B in RCGraph.all_rc_graphs(w, h):
                if len(B[-1]) != 0:
                    continue
                nxt = B.little_bump_desc()
                if nxt == B:
                    continue
                W, Wb = B.perm_word, nxt.perm_word
                if len(W) != len(Wb):
                    stats["length changed"] += 1
                    continue
                dec = [i for i in range(len(W)) if W[i] != Wb[i]]
                stats[f"letters decremented = {len(dec)}"] += 1

                rl1, rl2 = rootlist(W), rootlist(Wb)
                same = rl1 == rl2
                stats["rootlist unchanged" if same else "rootlist changed"] += 1

                # do the bare elements survive?
                bare1 = {r for r in rl1 if r[1] == 0}
                bare2 = {r for r in rl2 if r[1] == 0}
                stats["bare part unchanged" if bare1 == bare2 else "bare part changed"] += 1

                # are the changes confined to the values of decremented letters?
                touched_vals = {W[i] for i in dec} | {Wb[i] for i in dec}
                sym = (set(rl1) ^ set(rl2))
                if all(r[0] in touched_vals for r in sym):
                    stats["changes confined to decremented values"] += 1
                else:
                    stats["changes ESCAPE decremented values"] += 1
                    if len(examples) < 6:
                        examples.append((W, Wb, dec, window(rl1), window(rl2)))

    print(f"n = {n}  (single bumps)")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if examples:
        print("\nexamples where the rootlist change escapes the decremented values:")
        for W, Wb, dec, r1, r2 in examples:
            print(f"  {W} -> {Wb}  (positions {dec})")
            print(f"    rl before: {r1}")
            print(f"    rl after : {r2}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

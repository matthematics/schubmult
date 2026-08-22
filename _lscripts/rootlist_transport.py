"""Rootlist simulation via Nadeau-Tewari Lemma 5.7, and the base-case test.

rl(empty) is the set of bare integers; inserting a letter a deletes the two
rootlist elements bracketing a and inserts a. Elements of Z u Z are encoded as
pairs (value, tag) in lexicographic order, bare integers being (i, 0) and the
j-th occurrence of the value i in an injectified word being (i, j).

First the implementation is validated against forest equivalence: by Proposition
5.8, an adjacent transposition at p is a forest equivalence exactly when the two
letters are separated by at least two elements of the rootlist of the prefix.

Then, on the base-case instances (those where the bump decrements nothing before
position p, so that the upstairs and downstairs prefixes coincide), we ask
whether the same separating pair continues to separate the shifted letters.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph

LO, HI = -4, 40


def injectify(word):
    """Tag the j-th occurrence of a value with j, starting at 1."""
    seen = Counter()
    out = []
    for v in word:
        seen[v] += 1
        out.append((v, seen[v]))
    return out


def rootlist(prefix_letters):
    """rl of the prefix, by Lemma 5.7, restricted to a finite window."""
    rl = sorted((i, 0) for i in range(LO, HI))
    for a in prefix_letters:
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        if not below or not above:
            raise RuntimeError("window too small")
        r1, r2 = below[-1], above[0]
        rl = sorted(set(rl) - {r1, r2} | {a})
    return rl


def separators(a, b, rl):
    """The rootlist elements strictly between a and b."""
    lo, hi = (a, b) if a < b else (b, a)
    return [r for r in rl if lo < r < hi]


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()
    mismatches = []
    broken = []

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

                inj = injectify(d1)
                rlU = rootlist(inj[:p])
                seps = separators(inj[p], inj[p + 1], rlU)
                predicted = len(seps) >= 2
                actual = R1.forest_invariant == R2.forest_invariant
                stats[("validation", "agree" if predicted == actual else "MISMATCH")] += 1
                if predicted != actual and len(mismatches) < 5:
                    mismatches.append((d1, d2, p, len(seps), actual))

                if not actual:
                    continue

                # base case: the bump leaves the prefix alone
                up = B1.perm_word
                if any(up[i] != d1[i] for i in range(p)):
                    stats["feq: prefix touched"] += 1
                    continue
                stats["feq: base case (prefix untouched)"] += 1

                inj_up = injectify(up)
                rlU_up = rootlist(inj_up[:p])
                if rlU_up != rlU:
                    stats["  base case: rootlists differ anyway"] += 1
                seps_up = separators(inj_up[p], inj_up[p + 1], rlU_up)
                stats["  base case: upstairs separated" if len(seps_up) >= 2
                      else "  base case: upstairs NOT separated"] += 1
                # do the same two separators survive?
                survive = [r for r in seps if r in seps_up]
                stats[f"  base case: original separators surviving = {len(survive)}"] += 1
                if len(seps_up) < 2 and len(broken) < 5:
                    broken.append((up, d1, p, seps, seps_up))

    print(f"n = {n}")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if mismatches:
        print("\nvalidation mismatches:")
        for d1, d2, p, ns, actual in mismatches:
            print(f"  {d1} vs {d2} at p={p}: separators={ns}, forest-equivalent={actual}")
    if broken:
        print("\nbase-case failures:")
        for up, d1, p, s, su in broken:
            print(f"  up={up} down={d1} p={p}")
            print(f"    downstairs separators={s}  upstairs separators={su}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

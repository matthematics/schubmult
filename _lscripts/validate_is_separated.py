"""Validate LBS.is_separated against forest equivalence.

The prefix U of Nadeau-Tewari Proposition 5.8 is, in the reading order of the RC
graph, the reversal of the letters lying after the interchanged pair; letters are
tagged in insertion order, matching RCGraph.omega_invariant.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.combinatorics.rc_graph import RCGraph


def tagged_reversed(word):
    counts = {}
    out = []
    for a in reversed(word):
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def main(n):
    stats = Counter()
    bad = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = max(w.descents()) + 1
        for h in range(md, n + 1):
            rcs = list(RCGraph.all_rc_graphs(w, h))
            for R1, R2 in itertools.permutations(rcs, 2):
                d1, d2 = R1.perm_word, R2.perm_word
                diff = [i for i in range(len(d1)) if d1[i] != d2[i]]
                if len(diff) != 2 or diff[1] != diff[0] + 1:
                    continue
                p = diff[0]
                if not (d1[p] == d2[p + 1] and d1[p + 1] == d2[p]):
                    continue

                rev = tagged_reversed(d1)
                k = len(d1)
                # in reversed order the pair sits at the end of the prefix
                prefix = rev[: k - p - 2]
                b_letter = rev[k - p - 2]
                a_letter = rev[k - p - 1]

                res = omega_insertion(prefix)
                if res is None:
                    stats["insertion failed"] += 1
                    continue
                P = res[0]
                predicted = P.is_separated(a_letter, b_letter)
                actual = R1.forest_invariant == R2.forest_invariant
                stats["agree" if predicted == actual else "MISMATCH"] += 1
                if predicted != actual and len(bad) < 5:
                    bad.append((d1, d2, p, P.separators(a_letter, b_letter), actual))

    print(f"n = {n}")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if bad:
        print("\nmismatches:")
        for d1, d2, p, s, actual in bad:
            print(f"  {d1} vs {d2} at p={p}: separators={s} forest-equivalent={actual}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

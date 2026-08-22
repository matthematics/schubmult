"""Test: commutation-equivalent words have equal permutations after Z.

Claim (permutation form): if B1, B2 are bounded RC graphs of the same height whose
words lie in the same commutation class, then wof(Z(B1)) = wof(Z(B2)).

Commutation classes are computed by the recursive adjacent-commuting-swap closure.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def commutation_class(word):
    """Closure of `word` under swaps of adjacent letters differing by >= 2."""
    seen = {tuple(word)}
    stack = [tuple(word)]
    while stack:
        u = stack.pop()
        for i in range(len(u) - 1):
            if abs(u[i] - u[i + 1]) >= 2:
                v = u[:i] + (u[i + 1], u[i]) + u[i + 2:]
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
    return seen


def main(n):
    stats = Counter()
    failures = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        maxd = max(w.descents()) + 1
        for h in range(maxd, n + 1):
            graphs = [B for B in RCGraph.all_rc_graphs(w, h) if len(B[-1]) == 0]
            images = {}
            for B in graphs:
                try:
                    images[B] = B.zero_out_last_row()
                except Exception:
                    pass
            # group graphs by commutation class of their word
            classes = {}
            for B in images:
                cc = frozenset(commutation_class(B.perm_word))
                classes.setdefault(cc, []).append(B)
            for cc, members in classes.items():
                if len(members) < 2:
                    continue
                perms = {images[B].perm for B in members}
                trivial = maxd <= h - 1
                tag = "trivial" if trivial else "bump"
                stats[(tag, "pairs")] += len(members) * (len(members) - 1) // 2
                if len(perms) == 1:
                    stats[(tag, "same_perm_class")] += 1
                else:
                    stats[(tag, "DIFFERENT_perms")] += 1
                    if not trivial and len(failures) < 5:
                        failures.append((w, h, [B.perm_word for B in members],
                                         [images[B].perm for B in members]))

    print(f"n = {n}")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if failures:
        print("\nsample commutation classes with differing image permutations:")
        for w, h, words, perms in failures:
            print(f"  w={w} h={h}")
            print(f"    words={words}")
            print(f"    image perms={perms}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

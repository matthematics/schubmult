"""Probe the structure of forest (Omega) equivalence on reduced words."""

import sys
from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion


def pair_label(word):
    counts = {}
    out = []
    for a in word:
        a = int(a)
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def finv(word):
    return omega_insertion(pair_label(tuple(reversed(word))))[0]


def all_reduced_words(w):
    if w.inv == 0:
        return [()]
    out = []
    n = len(w)
    for i in range(1, n):
        if w[i - 1] > w[i]:
            v = w.swap(i - 1, i)
            out.extend([(i, *rest) for rest in all_reduced_words(v)])
    return out


def main(N=5):
    q3_fail = 0
    q3_tot = 0
    q1_fail = 0
    q1_tot = 0
    examples = []
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            words = [tuple(x) for x in all_reduced_words(w)]
            if not words:
                continue
            cls = defaultdict(list)
            for x in words:
                cls[finv(x)].append(x)
            for x in words:
                for i in range(len(x) - 1):
                    if abs(x[i] - x[i + 1]) >= 2:
                        y = x[:i] + (x[i + 1], x[i]) + x[i + 2 :]
                        q3_tot += 1
                        if finv(y) != finv(x):
                            q3_fail += 1
                            if len(examples) < 6:
                                examples.append(("commchange", tuple(w[:n]), x, y))
            for key, blk in cls.items():
                q1_tot += 1
                blkset = set(blk)
                seen = {blk[0]}
                stack = [blk[0]]
                while stack:
                    x = stack.pop()
                    for i in range(len(x) - 1):
                        if abs(x[i] - x[i + 1]) >= 2:
                            y = x[:i] + (x[i + 1], x[i]) + x[i + 2 :]
                            if y in blkset and y not in seen:
                                seen.add(y)
                                stack.append(y)
                if seen != blkset:
                    q1_fail += 1
                    if len(examples) < 12:
                        examples.append(("disconnected", tuple(w[:n]), sorted(blkset)))
    print(f"N={N}: commutation-preserves-forest: {q3_tot - q3_fail}/{q3_tot} (failures {q3_fail})")
    print(f"       classes connected by internal commutations: {q1_tot - q1_fail}/{q1_tot}")
    for e in examples:
        print("  ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

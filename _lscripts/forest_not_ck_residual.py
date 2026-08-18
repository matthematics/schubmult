"""Characterize the commutations that preserve the forest class but are NOT
Coxeter-Knuth moves (the residual left uncovered by Hamaker-Young).

Hypothesis to test: they are exactly the "isolated" commutations, where the two
swapped letters a, b have no letter strictly between them in value anywhere in
the word -- so no CK context exists.
"""

import sys
from collections import Counter
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.combinatorics.nilplactic import NilPlactic


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


def egP(word):
    return NilPlactic().ed_insert(*word)


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


def main(N=6):
    rows = []
    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            if w in seen:
                continue
            seen.add(w)
            for x in all_reduced_words(w):
                x = tuple(x)
                for i in range(len(x) - 1):
                    if abs(x[i] - x[i + 1]) < 2:
                        continue
                    y = x[:i] + (x[i + 1], x[i]) + x[i + 2 :]
                    if finv(x) != finv(y):
                        continue
                    if egP(x) == egP(y):
                        continue
                    rows.append((tuple(w[:n]), x, y, i))
    print(f"N={N}: forest-preserving but NOT Coxeter-Knuth: {len(rows)} instances")

    # hypothesis 1: no letter of the word lies strictly between the two swapped values
    h1 = 0
    # hypothesis 2: the word has length 2
    h2 = 0
    # hypothesis 3: every pair of letters in the word commutes (fully commutative word)
    h3 = 0
    lens = Counter()
    viol = []
    for wperm, x, y, i in rows:
        a, b = x[i], x[i + 1]
        lo, hi = min(a, b), max(a, b)
        lens[len(x)] += 1
        between = any(lo < c < hi for c in x)
        if not between:
            h1 += 1
        if len(x) == 2:
            h2 += 1
        if all(abs(p - q) >= 2 for j, p in enumerate(x) for q in x[j + 1 :]):
            h3 += 1
        else:
            if len(viol) < 8:
                viol.append((wperm, x, y))
    print(f"  word lengths: {dict(sorted(lens.items()))}")
    print(f"  no letter strictly between the swapped values : {h1}/{len(rows)}")
    print(f"  word has length exactly 2                      : {h2}/{len(rows)}")
    print(f"  word is fully commutative (all letters commute): {h3}/{len(rows)}")
    for v in viol:
        print("   not fully commutative:", v)
    print()
    print("  sample instances:")
    for r in rows[:10]:
        print("   ", r)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

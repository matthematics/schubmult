"""Is "commutation preserves the forest class" the same as "commutation is a
Coxeter-Knuth move"?

Hamaker-Young (Theorem littlecrystal in the paper) says a Little bump preserves
the Edelman-Greene recording tableau Q.  Coxeter-Knuth equivalence of reduced
words is the same as having a common EG insertion tableau P.  A BARE adjacent
commutation (e.g. 13 -> 31) is NOT in general a CK move; a commutation IN
CONTEXT (e.g. 132 -> 312) can be.

If the commutations that preserve the FOREST class are all CK moves, then
Hamaker-Young applies to exactly the moves the Gamma-ideal argument needs.

We cross-tabulate, for every adjacent commutation in every reduced word:
    forest-preserving?   vs   CK (EG P-tableau preserving)?
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
    """forest P-symbol, using the paper's convention (insert the reversed word)"""
    return omega_insertion(pair_label(tuple(reversed(word))))[0]


def egP(word):
    """Edelman-Greene insertion tableau"""
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


def main(N=5):
    tab = Counter()
    ex_forest_not_ck = []
    ex_ck_not_forest = []
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
                    f = finv(x) == finv(y)
                    c = egP(x) == egP(y)
                    tab[(f, c)] += 1
                    if f and not c and len(ex_forest_not_ck) < 6:
                        ex_forest_not_ck.append((tuple(w[:n]), x, y))
                    if c and not f and len(ex_ck_not_forest) < 6:
                        ex_ck_not_forest.append((tuple(w[:n]), x, y))
    print(f"N={N}  cross-tabulation over adjacent commutations")
    print(f"  forest-preserving AND CK      : {tab[(True, True)]}")
    print(f"  forest-preserving, NOT CK     : {tab[(True, False)]}   <-- breaks the Hamaker-Young route")
    print(f"  CK, NOT forest-preserving     : {tab[(False, True)]}")
    print(f"  neither                       : {tab[(False, False)]}")
    for e in ex_forest_not_ck:
        print("   forest but not CK:", e)
    for e in ex_ck_not_forest:
        print("   CK but not forest:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

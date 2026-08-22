"""Test whether the downstairs swap in eq:forestsimple_positions is a Coxeter-Knuth move.

By Edelman-Greene (Theorem 3.1 of Hamaker-Young), two reduced words are
Coxeter-Knuth equivalent iff they have the same EG insertion tableau P.
If forest equivalence of the two images forces P-equality, then Hamaker-Young
Lemma 3.3 (Little bumps commute with Coxeter-Knuth moves) proves the claim.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph
from schubmult.combinatorics.nilplactic import NilPlactic


def swap_positions(w1, w2):
    if len(w1) != len(w2) or w1 == w2:
        return None
    diff = [i for i in range(len(w1)) if w1[i] != w2[i]]
    if len(diff) != 2:
        return None
    p, q = diff
    if q != p + 1:
        return None
    if w1[p] == w2[q] and w1[q] == w2[p]:
        return p
    return None


def eg_p(word):
    return NilPlactic.from_word(tuple(word))


def ck_witness(word, p):
    """Is the swap at p,p+1 a Coxeter-Knuth move of type 1 or 2 (has a witness)?"""
    a, b = word[p], word[p + 1]
    lo, hi = min(a, b), max(a, b)
    # type 1: acb <-> cab, witness follows at p+2, strictly between
    if p + 2 < len(word) and lo < word[p + 2] < hi:
        return "type1"
    # type 2: bac <-> bca, witness precedes at p-1, strictly between
    if p - 1 >= 0 and lo < word[p - 1] < hi:
        return "type2"
    return None


def main(n):
    stats = Counter()
    mismatch = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        maxd = max(w.descents()) + 1
        for h in range(maxd, n + 1):
            if maxd <= h - 1:
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
                p = swap_positions(R1.perm_word, R2.perm_word)
                if p is None:
                    continue
                feq = R1.forest_invariant == R2.forest_invariant
                ok = swap_positions(B1.perm_word, B2.perm_word) == p
                same_P = eg_p(R1.perm_word) == eg_p(R2.perm_word)
                wit = ck_witness(R1.perm_word, p)
                key = ("feq" if feq else "noteq", "ok" if ok else "FAIL",
                       "sameP" if same_P else "diffP", wit or "no_witness")
                stats[key] += 1
                if feq and not same_P and len(mismatch) < 3:
                    mismatch.append((R1.perm_word, R2.perm_word, p))

    print(f"n = {n}  (nontrivial bump cases)")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if mismatch:
        print("\nforest-equivalent but different EG P-tableau:")
        for w1, w2, p in mismatch:
            print(f"  {w1} vs {w2} at p={p}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

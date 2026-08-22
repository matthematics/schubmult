"""Split the commutation moves by whether they are Coxeter-Knuth moves.

A commutation ab <-> ba is part of a Coxeter-Knuth move of type 1 (acb <-> cab)
or type 2 (bac <-> bca) when an adjacent third letter lies strictly between a and
b. Hamaker-Young Lemma 3.3 says Coxeter-Knuth moves commute with Little bumps, so
for those the terminal permutations should agree automatically; the question is
whether the bare commutations are exactly the ones needing the hypothesis.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def comm_at(u, v):
    if len(u) != len(v) or u == v:
        return None
    diff = [i for i in range(len(u)) if u[i] != v[i]]
    if len(diff) != 2:
        return None
    p, q = diff
    if q == p + 1 and u[p] == v[q] and u[q] == v[p] and abs(u[p] - u[q]) >= 2:
        return p
    return None


def ck_witness(word, p):
    a, b = word[p], word[p + 1]
    lo, hi = min(a, b), max(a, b)
    if p + 2 < len(word) and lo < word[p + 2] < hi:
        return "type1"
    if p - 1 >= 0 and lo < word[p - 1] < hi:
        return "type2"
    return "bare"


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            if md <= h - 1:
                continue
            reps = {}
            for B in RCGraph.all_rc_graphs(w, h):
                if len(B[-1]) != 0:
                    continue
                reps.setdefault(B.perm_word, B)
            imgs = {}
            for word, B in reps.items():
                try:
                    imgs[word] = B.zero_out_last_row().perm
                except Exception:
                    pass
            for u, v in itertools.combinations(imgs, 2):
                p = comm_at(u, v)
                if p is None:
                    continue
                kind = ck_witness(u, p)
                agree = imgs[u] == imgs[v]
                stats[(kind, "endpoints agree" if agree else "endpoints DIFFER")] += 1

    print(f"n = {n}  (nontrivial bump case, commutation moves)")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

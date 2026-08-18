"""Little bumps at the last-descent crossing, and the forest invariant.

Definitions follow writing/schubert_bialgebra_v2.tex:
  - word(R): concatenation of the rows of R (letters i+j-1).
  - lmap(R): single Little bump at the crossing whose root is (n,n+1), where
    n = maxd(w_R).  On words this is little_p(word(R)) with p the unique
    position such that deleting a_p leaves a reduced word for w s_n
    (strong exchange).
  - Z(R): iterate lmap until maxd <= n-1.

Checks:
  (A) word(Z(R)) equals the iterated bump of word(R)                [sanity]
  (B) single-bump lift: for B,B' in RC_h(w)^0 with maxd(w)=h,
          finv(lmap B) = finv(lmap B')  ==>  finv(B) = finv(B')
  (C) the number of bumps performed by Z is constant on RC_h(w)^0   [?]
"""

import sys
from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.indexed_forests import letterpair, omega_insertion
from schubmult.combinatorics.rc_graph import RCGraph


def pair_label(word):
    counts = {}
    out = []
    for a in word:
        a = int(a)
        counts[a] = counts.get(a, 0) + 1
        out.append(letterpair(a, counts[a]))
    return tuple(out)


def finv_word(word):
    return omega_insertion(pair_label(tuple(reversed(word))))[0]


def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])


def is_reduced(word):
    return perm_of(word).inv == len(word)


def maxd(perm):
    return len(perm.trimcode)


def little(word, p):
    """Little bump at position p (0-indexed)."""
    word = list(word)
    while True:
        if word[p] > 1:
            word[p] -= 1
        else:
            word = [a + 1 if q != p else a for q, a in enumerate(word)]
        if is_reduced(word):
            return tuple(word)
        # unique j != p whose deletion leaves a reduced word
        cands = [j for j in range(len(word)) if j != p and is_reduced(word[:j] + word[j + 1 :])]
        assert len(cands) == 1, (word, p, cands)
        p = cands[0]


def bump_position(word):
    """Position of the crossing with root (n,n+1), n = maxd."""
    w = perm_of(word)
    n = maxd(w)
    target = w * Permutation.ref_product(n)
    cands = [j for j in range(len(word)) if perm_of(word[:j] + word[j + 1 :]) == target and is_reduced(word[:j] + word[j + 1 :])]
    assert len(cands) == 1, (word, n, cands)
    return cands[0]


def bump(word):
    return little(word, bump_position(word))


def zero_word(word, h):
    """word(Z(R)) for R of height h with empty last row: bump until maxd <= h-1."""
    steps = 0
    while maxd(perm_of(word)) > h - 1:
        word = bump(word)
        steps += 1
    return word, steps


def rc0(perm, h):
    out = []
    for R in RCGraph.all_rc_graphs(perm, h):
        if len(R) == h and len(R[-1]) == 0:
            out.append(R)
    return out


def main(N=5):
    sanity_fail = 0
    sanity_tot = 0
    lift1_fail = 0
    lift1_tot = 0
    steps_vary = 0
    examples = []
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            d = maxd(w)
            if d < 1:
                continue
            for h in range(d, n + 1):
                graphs = rc0(w, h)
                if not graphs:
                    continue
                stepset = set()
                bumped = {}
                for B in graphs:
                    word = B.perm_word
                    zw, steps = zero_word(word, h)
                    stepset.add(steps)
                    sanity_tot += 1
                    if tuple(B.zero_out_last_row().perm_word) != tuple(zw):
                        sanity_fail += 1
                        if len(examples) < 5:
                            examples.append(("sanity", tuple(w[:n]), h, word, zw, tuple(B.zero_out_last_row().perm_word)))
                    if d == h:
                        bumped[B] = bump(word)
                if len(stepset) > 1:
                    steps_vary += 1
                if d == h and bumped:
                    grp = defaultdict(set)
                    for B, bw in bumped.items():
                        grp[finv_word(bw)].add(finv_word(B.perm_word))
                    for k, v in grp.items():
                        lift1_tot += 1
                        if len(v) > 1:
                            lift1_fail += 1
    print(f"N={N}: sanity word(Z)=iterated bump: {sanity_tot - sanity_fail}/{sanity_tot}")
    print(f"       single-bump lift: {lift1_tot - lift1_fail}/{lift1_tot} classes (failures {lift1_fail})")
    print(f"       (w,h) with non-constant bump count: {steps_vary}")
    for e in examples:
        print("  ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

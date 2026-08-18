"""Step (2), tested on ALL single-crossing anti-diagonal slides.

Earlier I verified the intertwining only for pairs whose WORDS are one adjacent
commutation apart.  Step (1) (connectivity of forest classes) uses slides of
arbitrary displacement, which need not be word-commutations.  So the two steps
were quantified over different sets.  This script closes that mismatch.

A SLIDE: R1 -> R2, both RC graphs of the same permutation and height, with
exactly one crossing moved, from (i,j) to (i',j'), with i+j = i'+j' (so the
letter L = i+j-1 is preserved).

Tested:
  (2a) if R1, R2 both lie in the image of Z, the fibers are indexed by the same
       permutation set;
  (2b) for each w, the matched pair B_w, B'_w differs by the IDENTICAL slide
       (same letter, same source row, same target row);
  (2c) how many of these slides are word-commutations (the set tested before).
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc_all(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h]


def rc0(perm, h):
    return [R for R in rc_all(perm, h) if len(R[-1]) == 0]


def crossings(R):
    out = set()
    for i, row in enumerate(R, start=1):
        for L in row:
            out.add((i, L - i + 1))
    return out


def slide_sig(R1, R2):
    """(letter, i, i') if R1->R2 is a single anti-diagonal slide, else None."""
    c1, c2 = crossings(R1), crossings(R2)
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    if i + j != ip + jp:
        return None
    return (i + j - 1, i, ip)


def comm_neighbors(word):
    out = []
    for i in range(len(word) - 1):
        if abs(word[i] - word[i + 1]) >= 2:
            out.append(word[:i] + (word[i + 1], word[i]) + word[i + 2 :])
    return out


def main(N=5):
    idx_tot = idx_bad = 0
    lift_tot = lift_bad = 0
    disp_hist = Counter()
    iscomm = Counter()
    bad_examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            u = Permutation(list(pw))
            if maxd(u) == 0:
                continue
            for m in range(max(maxd(u), 1), n + 1):
                if (u, m) in seen:
                    continue
                seen.add((u, m))

                fiber = defaultdict(dict)
                for B in rc0(u, m + 1):
                    fiber[B.zero_out_last_row()][B.perm] = B

                Rs = list(fiber)
                for R1 in Rs:
                    for R2 in Rs:
                        if R1 is R2:
                            continue
                        sig = slide_sig(R1, R2)
                        if sig is None:
                            continue
                        L, i, ip = sig
                        disp_hist[ip - i] += 1
                        iscomm[tuple(R2.perm_word) in comm_neighbors(tuple(R1.perm_word))] += 1

                        idx_tot += 1
                        if set(fiber[R1]) != set(fiber[R2]):
                            idx_bad += 1
                            continue
                        for w in fiber[R1]:
                            B, Bp = fiber[R1][w], fiber[R2][w]
                            lift_tot += 1
                            su = slide_sig(B, Bp)
                            if su != sig:
                                lift_bad += 1
                                if len(bad_examples) < 6:
                                    bad_examples.append(
                                        (
                                            tuple(u[:n]),
                                            m,
                                            tuple(w[: n + 1]),
                                            f"down {sig}",
                                            f"up {su}",
                                            [tuple(r) for r in R1],
                                            [tuple(r) for r in R2],
                                            [tuple(r) for r in B],
                                            [tuple(r) for r in Bp],
                                        )
                                    )
    print(f"N={N}: slides between graphs in the image of Z = {idx_tot}")
    print(f"  displacement histogram = {dict(sorted(disp_hist.items()))}")
    print(f"  of these, word-commutations: {iscomm[True]}   NOT word-commutations: {iscomm[False]}")
    print(f"  fibers indexed by same permutation set: {idx_tot - idx_bad}/{idx_tot} (bad {idx_bad})")
    print(f"  (2b) upstairs is the IDENTICAL slide: {lift_tot - lift_bad}/{lift_tot} (bad {lift_bad})")
    for e in bad_examples:
        print("   FAIL:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

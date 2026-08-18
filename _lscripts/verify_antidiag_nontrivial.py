"""CORRECTED test of Lemma zeroantidiag, on the NONTRIVIAL branch of Algorithm 3.

BUG IN THE EARLIER SCRIPTS.  They looped
        for m in range(max(maxd(u), 1), n + 1)
and built the fiber from rc0(u, m + 1), i.e. graphs of height m+1 with
permutation u.  That forces maxd(u) <= m, which is precisely the condition
under which the FIRST BRANCH of Algorithm 3 applies: Z merely deletes the empty
last row.  So B was always R with a row appended, and the intertwining was true
by construction.  Every earlier count (42188/42188 etc.) was vacuous.

Here we take h = maxd(w) exactly, the case in which Algorithm 3 actually runs
its insertion loop.  Upstairs: B in RC_h(w)^0 with maxd(w) = h.  Downstairs:
R = Z(B) of height m = h - 1.

We test, for pairs R1, R2 of the SAME permutation differing by one antidiagonal
move sigma = (letter L, row i -> row i'):
  (a) is sigma applicable to B_1 at all (crossing present, target free)?
  (b) is B_2 = sigma(B_1)?
  (c) if not, what IS the relation between B_1 and B_2 -- same rows with a
      different letter (the user's sigma'), different rows, or not a single
      antidiagonal move at all?
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def crossings(R):
    out = set()
    for i, row in enumerate(R, start=1):
        for L in row:
            out.add((i, L - i + 1))
    return out


def slide_sig(c1, c2):
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    if i + j != ip + jp:
        return None
    return (i + j - 1, i, ip)


def main(N=6):
    pairs = 0
    trivial_guard = 0
    applicable = 0
    equals_sigma = 0
    kind = Counter()
    letter_shift = Counter()
    rowpat = Counter()
    examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2:
                continue
            if (w, h) in seen:
                continue
            seen.add((w, h))

            graphs = rc0(w, h)
            if len(graphs) < 2:
                continue

            Zinv = {}
            for B in graphs:
                R = B.zero_out_last_row()
                # sanity: confirm Z is NOT the trivial row deletion here
                if crossings(R) == crossings(B):
                    trivial_guard += 1
                Zinv[R] = B

            Rs = list(Zinv)
            for R1 in Rs:
                c1 = crossings(R1)
                for R2 in Rs:
                    if R1 is R2:
                        continue
                    if R1.perm != R2.perm:
                        continue
                    sig = slide_sig(c1, crossings(R2))
                    if sig is None:
                        continue
                    L, i, ip = sig
                    pairs += 1

                    B1, B2 = Zinv[R1], Zinv[R2]
                    d1, d2 = crossings(B1), crossings(B2)
                    j, jp = L - i + 1, L - ip + 1

                    if ((i, j) in d1) and ((ip, jp) not in d1):
                        applicable += 1
                        if (d1 - {(i, j)}) | {(ip, jp)} == d2:
                            equals_sigma += 1

                    su = slide_sig(d1, d2)
                    if su is None:
                        kind["not a single antidiagonal move"] += 1
                        if len(examples) < 6:
                            examples.append(
                                (tuple(w[:n]), h, f"down {sig}", "up: multi-crossing",
                                 [tuple(r) for r in R1], [tuple(r) for r in R2],
                                 [tuple(r) for r in B1], [tuple(r) for r in B2])
                            )
                    else:
                        Lu, iu, ipu = su
                        letter_shift[Lu - L] += 1
                        rowpat[(iu - i, ipu - ip)] += 1
                        if su == sig:
                            kind["identical move"] += 1
                        elif (iu, ipu) == (i, ip):
                            kind["same rows, different letter"] += 1
                            if len(examples) < 6:
                                examples.append(
                                    (tuple(w[:n]), h, f"down {sig}", f"up {su}",
                                     [tuple(r) for r in R1], [tuple(r) for r in R2],
                                     [tuple(r) for r in B1], [tuple(r) for r in B2])
                                )
                        else:
                            kind["different rows"] += 1
                            if len(examples) < 6:
                                examples.append(
                                    (tuple(w[:n]), h, f"down {sig}", f"up {su}",
                                     [tuple(r) for r in R1], [tuple(r) for r in R2],
                                     [tuple(r) for r in B1], [tuple(r) for r in B2])
                                )

    print(f"N={N}  (h = maxd(w), the NONTRIVIAL branch of Algorithm 3)")
    print(f"  instances where Z was still the trivial row deletion: {trivial_guard}")
    print(f"  same-permutation antidiagonal pairs downstairs: {pairs}")
    print()
    print(f"  (a) sigma applicable to B_1: {applicable}/{pairs}")
    print(f"  (b) B_2 == sigma(B_1):       {equals_sigma}/{pairs}")
    print(f"  (c) upstairs relation: {dict(kind)}")
    print(f"      letter shift L_up - L_down: {dict(sorted(letter_shift.items()))}")
    print(f"      (i_up - i, i'_up - i'):     {dict(sorted(rowpat.items()))}")
    for e in examples:
        print("\n   MISMATCH:")
        print(f"     w={e[0]} h={e[1]}  {e[2]}   {e[3]}")
        print(f"     R1={e[4]}  R2={e[5]}")
        print(f"     B1={e[6]}  B2={e[7]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

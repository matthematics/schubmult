"""Strengthen step (2) to n=6, and test the proof strategy.

PROOF STRATEGY for step (2).  Let sigma be a slide, B in RC_{m+1}(w)^0 with
Z(B) = R1, and R2 = sigma(R1).  A slide preserves the permutation, so
w_{sigma(B)} = w.  Suppose we know

        trm( Z(sigma(B)) )  =  trm( R2 ).                     (*)

Then Z(sigma(B)) and R2 both have height m, the same trim, and the same weight
(Z preserves weight by lemma:zeroweight, and both sides have weight
wt(R1) - e_i + e_{i'}).  By theorem:modactrc + lemma:factor both satisfy
w' down_1 u' with l(w') - l(u') = p, and both are down_{m+1}-below w; so both
lie in the set W of lemma:transition_interchange, which has at most one element
by lemma:transition_unique.  Hence they have the same permutation, and by
theorem:modactrc they are equal.

So everything reduces to (*).  By lemma:ztrim, trm(Z(X)) = Z(trm(X)), so
(*) reads   Z( trm(sigma(B)) ) = trm( sigma(R1) ).

  - If the slide avoids row 1, trm(sigma(X)) = sigma'(trm(X)) for the shifted
    slide sigma', and (*) is the SAME statement at height m: induction.
  - If the slide touches row 1, trm(sigma(X)) is trm(X) with a crossing ADDED
    or REMOVED, and (*) becomes an "insertion commutes with Z" statement.

This script tests the two reductions separately, so we know exactly which part
is genuinely open.
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
    c1, c2 = crossings(R1), crossings(R2)
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    if i + j != ip + jp:
        return None
    return (i + j - 1, i, ip)


def main(N=6):
    tot = bad = 0
    touch1_tot = touch1_bad = 0
    avoid1_tot = avoid1_bad = 0
    disp = Counter()
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
                        disp[ip - i] += 1
                        if set(fiber[R1]) != set(fiber[R2]):
                            continue
                        touches1 = (i == 1) or (ip == 1)
                        for w in fiber[R1]:
                            B, Bp = fiber[R1][w], fiber[R2][w]
                            tot += 1
                            ok = slide_sig(B, Bp) == sig
                            if touches1:
                                touch1_tot += 1
                                if not ok:
                                    touch1_bad += 1
                            else:
                                avoid1_tot += 1
                                if not ok:
                                    avoid1_bad += 1
                            if not ok:
                                bad += 1
                                if len(bad_examples) < 5:
                                    bad_examples.append((tuple(u[:n]), m, sig, slide_sig(B, Bp)))

    print(f"N={N}: slide instances = {tot}")
    print(f"  displacement histogram = {dict(sorted(disp.items()))}")
    print(f"  (2b) upstairs is the IDENTICAL slide: {tot - bad}/{tot} (bad {bad})")
    print(f"     slides AVOIDING row 1 (handled by the induction): {avoid1_tot - avoid1_bad}/{avoid1_tot}")
    print(f"     slides TOUCHING row 1 (the residual case):        {touch1_tot - touch1_bad}/{touch1_tot}")
    for e in bad_examples:
        print("   FAIL:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

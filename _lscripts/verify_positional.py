"""The positional approach to lemma:forestsimple.

Z preserves the weight vector (Lemma zeroweight), so word(B) and word(Z(B)) have
the same length AND the same row structure: position k of one corresponds to
position k of the other.  (Verified earlier: the row of every crossing is
preserved, 524/524 at n=5.)  Letters change, positions do not.

USER'S OBSERVATION.  Zeroing/Little bumping does not preserve commutation.  But
when it fails to, the PERMUTATION CHANGES: the two words are still distinct
outputs at the same two positions, but they are reduced words for different
permutations.

So the dichotomy at positions p, p+1 should be:
    letters downstairs still commute  =>  same permutation, still a commutation
    letters downstairs are adjacent   =>  different permutations

For lemma:forestsimple we need the BACKWARD form: R1, R2 have the SAME
permutation and differ by a swap at p,p+1; their preimages B1, B2 (common w)
should then differ by a swap at p,p+1, hence be forest equivalent.

Tests (nontrivial branch h = maxd(w) reported separately from the trivial one):
 (P) position correspondence: word(B) and word(Z(B)) have equal length and equal
     row-block structure
 (F) forward: word(B2) = word(B1) swapped at p,p+1  =>  word(Z(B2)) =
     word(Z(B1)) swapped at p,p+1
 (D) the dichotomy: letters downstairs commute <=> permutations agree
 (B) backward: R1, R2 same permutation, swap at p,p+1, preimages differ by a
     swap at p,p+1
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(p):
    return len(p.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def word(R):
    return tuple(R.perm_word)


def rowblocks(R):
    return tuple(len(r) for r in R)


def swap_at(w, p):
    return w[:p] + (w[p + 1], w[p]) + w[p + 2 :]


def swap_positions(w1, w2):
    """positions p with w2 = w1 swapped at p, p+1 (letters distinct)"""
    if len(w1) != len(w2):
        return []
    return [p for p in range(len(w1) - 1) if swap_at(w1, p) == w2]


def main(N=6):
    P_tot = P_ok = 0
    F = defaultdict(lambda: [0, 0])
    D = Counter()
    Bk = defaultdict(lambda: [0, 0])
    ex_F = []
    ex_B = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            for h in range(max(maxd(w), 1), n + 1):
                if (w, h) in seen:
                    continue
                seen.add((w, h))
                graphs = rc0(w, h)
                if len(graphs) < 2:
                    continue
                branch = "nontrivial" if maxd(w) == h else "trivial"

                Z = {}
                for B in graphs:
                    R = B.zero_out_last_row()
                    Z[B] = R
                    P_tot += 1
                    if len(word(B)) == len(word(R)) and rowblocks(B)[: h - 1] == rowblocks(R):
                        P_ok += 1

                # (F) and (D): forward
                for B1 in graphs:
                    for B2 in graphs:
                        if B1 is B2:
                            continue
                        ps = swap_positions(word(B1), word(B2))
                        if not ps:
                            continue
                        R1, R2 = Z[B1], Z[B2]
                        p = ps[0]
                        ok = swap_positions(word(R1), word(R2)) != []
                        F[branch][0] += 1
                        if ok:
                            F[branch][1] += 1
                        elif len(ex_F) < 4:
                            ex_F.append((branch, tuple(w[:n]), h, word(B1), word(B2), word(R1), word(R2)))
                        # dichotomy
                        if p + 1 < len(word(R1)):
                            a, b = word(R1)[p], word(R1)[p + 1]
                            commutes = abs(a - b) >= 2
                            sameperm = R1.perm == R2.perm
                            D[(branch, commutes, sameperm)] += 1

                # (B) backward
                byR = {}
                for B in graphs:
                    byR[Z[B]] = B
                Rs = list(byR)
                for R1 in Rs:
                    for R2 in Rs:
                        if R1 is R2 or R1.perm != R2.perm:
                            continue
                        if not swap_positions(word(R1), word(R2)):
                            continue
                        B1, B2 = byR[R1], byR[R2]
                        Bk[branch][0] += 1
                        if swap_positions(word(B1), word(B2)):
                            Bk[branch][1] += 1
                        elif len(ex_B) < 5:
                            ex_B.append((branch, tuple(w[:n]), h, word(R1), word(R2), word(B1), word(B2)))

    print(f"N={N}")
    print(f"  (P) position/row-block correspondence between B and Z(B): {P_ok}/{P_tot}")
    print()
    print("  (F) forward: swap upstairs => swap downstairs")
    for k, (t, o) in F.items():
        print(f"      {k}: {o}/{t}")
    print()
    print("  (D) dichotomy (branch, letters commute downstairs, permutations agree) -> count")
    for k in sorted(D, key=str):
        print(f"      {k}: {D[k]}")
    print()
    print("  (B) BACKWARD (what forestsimple needs): same perm + swap downstairs => swap upstairs")
    for k, (t, o) in Bk.items():
        print(f"      {k}: {o}/{t}")
    for e in ex_F:
        print("\n   FORWARD FAILS:", e)
    for e in ex_B:
        print("\n   BACKWARD FAILS:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

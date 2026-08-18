"""Test the v1 right-ideal argument, which never needs interior moves.

v1 (writing/schubert_bialgebra.tex, proof of lemma:forestideal):

  reduce to the case where a and b are the LAST TWO letters of word(R1), with b
  the last descent of w_{R1}, and R the empty graph.  Then:
    - "for any term not identical to R1, the suffix of the word will be a b'
       with b' > b";
    - "for each b' there will be precisely one matching term with suffix b' a
       and identical prefix";
    - "the local swap of a and b' is a valid local swap".

Taking R = the empty graph of height 1, Definition ringproduct + lemma:clipempty
give  R1 (+) empty_1 = sum_{Z(B) = R1} B, so the terms of the product ARE the
fibre of Z.

Claims tested, for pairs R1, R2 of the same permutation whose words differ by
swapping the last two letters:
  (C0) how often is such a pair forest equivalent
  (C1) every B in the fibre of R1 has word = word(R1) with only the LAST letter
       changed, and the new letter b' satisfies b' >= b
  (C2) the fibres of R1 and R2 are matched by b' <-> b', with words
       prefix.a.b'  <->  prefix.b'.a
  (C3) matched terms are forest equivalent
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(p):
    return len(p.trimcode)


def word(R):
    return tuple(R.perm_word)


def all_empty_last(h, n):
    """all bounded RC graphs of height h, empty last row, perms in S_n"""
    out = []
    for pw in permutations(range(1, n + 1)):
        w = Permutation(list(pw))
        if maxd(w) > h:
            continue
        for B in RCGraph.all_rc_graphs(w, h):
            if len(B) == h and len(B[-1]) == 0:
                out.append(B)
    return out


def main(N=5):
    c0_tot = c0_fe = 0
    c1_tot = c1_ok = 0
    c1_ge = 0
    c2_tot = c2_ok = 0
    c3_tot = c3_ok = 0
    lastletter_only = Counter()
    skipped_hyp = [0]
    examples = []

    for n in range(2, N + 1):
        for m in range(1, n + 1):
            graphs = all_empty_last(m + 1, n)
            if not graphs:
                continue
            fib = defaultdict(list)
            for B in graphs:
                fib[B.zero_out_last_row()].append(B)

            Rs = list(fib)
            byword = defaultdict(list)
            for R in Rs:
                byword[(R.perm, word(R))].append(R)

            for R1 in Rs:
                w1 = word(R1)
                if len(w1) < 2:
                    continue
                a, b = w1[-2], w1[-1]
                if abs(a - b) < 2:
                    continue
                # v1 hypothesis: b is the last descent of w_{R1}
                if b != maxd(R1.perm):
                    skipped_hyp[0] += 1
                    continue
                w2 = w1[:-2] + (b, a)
                for R2 in byword.get((R1.perm, w2), []):
                    c0_tot += 1
                    fe = R1.forest_invariant == R2.forest_invariant
                    if fe:
                        c0_fe += 1
                    if not fe:
                        continue

                    # (C1) fibre words
                    f1, f2 = fib[R1], fib[R2]
                    ok1 = True
                    for B in f1:
                        c1_tot += 1
                        wb = word(B)
                        if len(wb) == len(w1) and wb[:-1] == w1[:-1]:
                            c1_ok += 1
                            if wb[-1] >= b:
                                c1_ge += 1
                        else:
                            ok1 = False
                            if len(examples) < 4:
                                examples.append(("C1", tuple(R1.perm[:n]), m, w1, wb))
                    lastletter_only[ok1] += 1

                    # (C2) matching by the varying letter
                    key1 = {}
                    for B in f1:
                        wb = word(B)
                        if wb[:-1] == w1[:-1]:
                            key1[wb[-1]] = B
                    key2 = {}
                    for B in f2:
                        wb = word(B)
                        if len(wb) == len(w2) and wb[:-1] == w2[:-1]:
                            key2[wb[-1]] = B
                    # v1 predicts fibre of R2 has words prefix . b' . a
                    key2b = {}
                    for B in f2:
                        wb = word(B)
                        if len(wb) >= 2 and wb[:-2] == w1[:-2] and wb[-1] == a:
                            key2b[wb[-2]] = B
                    c2_tot += 1
                    if set(key1) == set(key2b) and len(key1) == len(f1) and len(key2b) == len(f2):
                        c2_ok += 1
                    # (C3) matched terms forest equivalent
                    for bp, B in key1.items():
                        if bp in key2b:
                            c3_tot += 1
                            if B.forest_invariant == key2b[bp].forest_invariant:
                                c3_ok += 1

    print(f"N={N}  (skipped, b not the last descent: {skipped_hyp[0]})")
    print(f"  (C0) pairs differing by a swap of the LAST TWO letters: {c0_tot}, forest equivalent: {c0_fe}")
    print()
    print(f"  (C1) fibre words change ONLY the last letter: {c1_ok}/{c1_tot}")
    print(f"       ... and the new letter is >= b:          {c1_ge}/{c1_tot}")
    print(f"       per pair, all terms conform: {dict(lastletter_only)}")
    print(f"  (C2) fibres matched by the varying letter, words prefix.a.b' <-> prefix.b'.a: {c2_ok}/{c2_tot}")
    print(f"  (C3) matched terms are forest equivalent: {c3_ok}/{c3_tot}")
    for e in examples:
        print("   ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

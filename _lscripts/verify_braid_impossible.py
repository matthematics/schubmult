"""The braid case cannot occur upstairs: a proof-shaped test.

Positions are preserved by Little bumps (verified: the transposed pair stays at
p, p+1, 654/654).  Suppose R1, R2 downstairs have the same permutation and words
differing by a transposition at p, p+1, and let B1, B2 be their preimages in
RC_h(w)^0 -- note BOTH have the SAME permutation w.

Then word(B1) and word(B2) differ by a transposition at p, p+1, carrying letters
a, b.  If |a - b| = 1 then s_a s_b != s_b s_a, so B1 and B2 would be reduced
words for DIFFERENT permutations, contradicting that both equal w.  Hence
|a - b| >= 2: the braid case is impossible upstairs, and the upstairs relation
is a commutation.

USER'S CLAIM, tested separately: braid-related words are never forest
equivalent.

Tests:
  (1) in every simple-case instance, do the upstairs words differ by a
      transposition at the SAME positions, with commuting letters?
  (2) braid-related reduced words: same permutation?  forest equivalent?
  (3) the remaining question: is the upstairs commutation always
      forest-permitted (i.e. are B1, B2 forest equivalent)?
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from verify_forest_littlebump import finv_word, is_reduced, maxd, perm_of


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def word(R):
    return tuple(R.perm_word)


def transpose_at(w, p):
    return w[:p] + (w[p + 1], w[p]) + w[p + 2 :]


def trans_positions(w1, w2):
    if len(w1) != len(w2):
        return []
    return [p for p in range(len(w1) - 1) if transpose_at(w1, p) == w2 and w1[p] != w1[p + 1]]


def main(N=6):
    simple = 0
    up_same_pos = 0
    up_commute = 0
    up_braid = 0
    up_forest = 0
    updown_pos = Counter()
    ex = []

    # user's claim
    braid_tot = 0
    braid_sameperm = 0
    braid_forest = 0

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2 or (w, h) in seen:
                continue
            seen.add((w, h))

            graphs = rc0(w, h)
            if len(graphs) < 2:
                continue
            Zinv = {}
            for B in graphs:
                Zinv[B.zero_out_last_row()] = B
            Rs = list(Zinv)

            # (2) braid-related reduced words for this permutation
            ws = sorted({word(B) for B in graphs})
            for U in ws:
                for p in range(len(U) - 1):
                    if abs(U[p] - U[p + 1]) != 1:
                        continue
                    V = transpose_at(U, p)
                    if not is_reduced(V):
                        continue
                    braid_tot += 1
                    if perm_of(V) == perm_of(U):
                        braid_sameperm += 1
                        if finv_word(U) == finv_word(V):
                            braid_forest += 1

            # (1),(3) the simple case
            for R1 in Rs:
                for R2 in Rs:
                    if R1 is R2 or R1.perm != R2.perm:
                        continue
                    ps = trans_positions(word(R1), word(R2))
                    if not ps:
                        continue
                    if R1.forest_invariant != R2.forest_invariant:
                        continue
                    simple += 1
                    B1, B2 = Zinv[R1], Zinv[R2]
                    u1, u2 = word(B1), word(B2)
                    qs = trans_positions(u1, u2)
                    if set(qs) & set(ps):
                        up_same_pos += 1
                    if qs:
                        q = qs[0]
                        updown_pos[tuple(sorted(set(qs) & set(ps))) != ()] += 1
                        if abs(u1[q] - u1[q + 1]) >= 2:
                            up_commute += 1
                        else:
                            up_braid += 1
                            if len(ex) < 5:
                                ex.append((tuple(w[:n]), h, word(R1), word(R2), u1, u2))
                    elif len(ex) < 5:
                        ex.append(("NOT A TRANSPOSITION", tuple(w[:n]), h, word(R1), word(R2), u1, u2))
                    if B1.forest_invariant == B2.forest_invariant:
                        up_forest += 1

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  simple-case instances: {simple}")
    print(f"    upstairs words differ by a transposition at the SAME position: {up_same_pos}/{simple}")
    print(f"    upstairs letters COMMUTE: {up_commute}   upstairs letters BRAID: {up_braid}")
    print(f"    upstairs forest equivalent (= lemma forestsimple): {up_forest}/{simple}")
    print()
    print("  USER'S CLAIM: braid-related reduced words")
    print(f"    braid transpositions that stay reduced: {braid_tot}")
    print(f"      ... with the SAME permutation: {braid_sameperm}")
    print(f"      ... and forest equivalent:      {braid_forest}")
    for e in ex:
        print("\n   ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

"""Classify adjacent-commutation moves on RC graphs into crystal moves and
"crystal breakers" (user's term), and test the lifting conjecture on each kind.

Conventions.  An RCGraph is a list of rows; row i is the tuple of LETTERS
(i+j-1) in that row, listed in decreasing order.  word(R) is the concatenation
of the rows.  A crystal operator e_i/f_i moves one crossing between rows i and
i+1, preserving the multiset of letters, hence permuting the word.

The user's example:
    A = {(1,1),(2,2)} as (row,column) crossings -> letters 1 (row 1), 3 (row 2)
        so A = RCGraph([(1,),(3,)]),  word 13, weight (1,1)
    B = row 1 carries letters 3 and 1
        so B = RCGraph([(3,1),()]),   word 31, weight (2,0)
"""

import sys
from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc_all(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h]


def rc0(perm, h):
    return [R for R in rc_all(perm, h) if len(R[-1]) == 0]


def comm_neighbors(word):
    out = []
    for i in range(len(word) - 1):
        if abs(word[i] - word[i + 1]) >= 2:
            out.append(word[:i] + (word[i + 1], word[i]) + word[i + 2 :])
    return out


def crystal_related(R1, R2):
    """Is R2 obtained from R1 by a single crystal operator e_i or f_i?"""
    for i in range(1, len(R1) + 1):
        for op in ("raising_operator", "lowering_operator"):
            try:
                got = getattr(R1, op)(i)
            except Exception:
                got = None
            if got is not None and got == R2:
                return (op, i)
    return None


def show_example():
    A = RCGraph([(1,), (3,)])
    B = RCGraph([(3, 1), ()])
    for name, R in (("A", A), ("B", B)):
        print(f"  {name}: rows={[tuple(r) for r in R]} word={tuple(R.perm_word)} wt={tuple(len(r) for r in R)} perm={tuple(R.perm[:4])}")
    print(f"  same permutation: {A.perm == B.perm}")
    print(f"  words one commutation apart: {tuple(B.perm_word) in comm_neighbors(tuple(A.perm_word))}")
    print(f"  crystal-related A->B: {crystal_related(A, B)}")
    print(f"  crystal-related B->A: {crystal_related(B, A)}")
    print(f"  same forest class: {A.forest_invariant == B.forest_invariant}")
    for i in (1,):
        print(f"  e_{i}(A)={A.raising_operator(i)}  f_{i}(A)={A.lowering_operator(i)}")


def main(N=5):
    print("=== user's example ===")
    show_example()

    print("\n=== classification of commutation moves ===")
    n_cry = n_break = 0
    cry_forest_ok = break_forest_ok = 0
    break_examples = []

    seen = set()
    # lifting test data
    lift = {"crystal": [0, 0], "breaker": [0, 0]}

    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            u = Permutation(list(pw))
            if maxd(u) == 0:
                continue
            for m in range(max(maxd(u), 1), n + 1):
                if (u, m) in seen:
                    continue
                seen.add((u, m))

                # fibers for the lifting test
                fiber = defaultdict(dict)
                for B in rc0(u, m + 1):
                    fiber[B.zero_out_last_row()][B.perm] = B

                graphs = rc_all(u, m)
                bywordset = {}
                for R in graphs:
                    bywordset.setdefault(tuple(R.perm_word), []).append(R)

                for R1 in graphs:
                    w1 = tuple(R1.perm_word)
                    for w2 in comm_neighbors(w1):
                        for R2 in bywordset.get(w2, []):
                            rel = crystal_related(R1, R2)
                            kind = "crystal" if rel else "breaker"
                            if rel:
                                n_cry += 1
                                if R1.forest_invariant == R2.forest_invariant:
                                    cry_forest_ok += 1
                            else:
                                n_break += 1
                                if R1.forest_invariant == R2.forest_invariant:
                                    break_forest_ok += 1
                                if len(break_examples) < 6:
                                    break_examples.append(
                                        (
                                            [tuple(r) for r in R1],
                                            [tuple(r) for r in R2],
                                            w1,
                                            w2,
                                            R1.forest_invariant == R2.forest_invariant,
                                        )
                                    )
                            # lifting test: only for graphs that are in the image of Z
                            if R1 in fiber and R2 in fiber and set(fiber[R1]) == set(fiber[R2]):
                                for w in fiber[R1]:
                                    b1, b2 = fiber[R1][w], fiber[R2][w]
                                    x, y = tuple(b1.perm_word), tuple(b2.perm_word)
                                    lift[kind][0] += 1
                                    if y in comm_neighbors(x):
                                        lift[kind][1] += 1

    print(f"  crystal moves: {n_cry}   of which forest-class-preserving: {cry_forest_ok}")
    print(f"  crystal breakers: {n_break}   of which forest-class-preserving: {break_forest_ok}")
    print("\n=== lifting conjecture, split by kind ===")
    for k, (tot, ok) in lift.items():
        print(f"  {k}: one commutation lifts to one commutation in {ok}/{tot}")
    print("\n=== sample crystal breakers (rows1, rows2, word1, word2, same forest class) ===")
    for e in break_examples:
        print("   ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

"""The bridging statement my architecture actually needs.

For a fiber-matched pair related by an antidiagonal move sigma:

    R2 = sigma(R1)   (height m)          B2 = sigma(B1)   (height m+1)

is it true that

    R1 =_forest R2   <==>   B1 =_forest B2   ?

The "<==" direction, combined with connectivity of forest classes by
antidiagonal moves, is what gives Lemma forestlift.  Chaining alone does NOT
give it: an antidiagonal move relocates a letter of the word, so it need not
preserve the forest class, and preservation downstairs does not formally imply
preservation upstairs.

We test both directions, and also record how often an antidiagonal move
preserves the forest class at all (to confirm the statement is not vacuous).
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


def main(N=5):
    tab = Counter()          # (down_same, up_same) -> count
    examples = defaultdict(list)

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
                        if set(fiber[R1]) != set(fiber[R2]):
                            continue
                        down = R1.forest_invariant == R2.forest_invariant
                        for w in fiber[R1]:
                            B1, B2 = fiber[R1][w], fiber[R2][w]
                            if slide_sig(B1, B2) != sig:
                                tab["MOVE-MISMATCH"] += 1
                                continue
                            up = B1.forest_invariant == B2.forest_invariant
                            tab[(down, up)] += 1
                            if down != up and len(examples[(down, up)]) < 4:
                                examples[(down, up)].append(
                                    (
                                        tuple(u[:n]),
                                        m,
                                        sig,
                                        [tuple(r) for r in R1],
                                        [tuple(r) for r in R2],
                                        [tuple(r) for r in B1],
                                        [tuple(r) for r in B2],
                                    )
                                )

    print(f"N={N}: (down preserves class, up preserves class) -> count")
    for k in [(True, True), (True, False), (False, True), (False, False), "MOVE-MISMATCH"]:
        if k in tab:
            print(f"   {k}: {tab[k]}")
    dt = tab[(True, True)] + tab[(True, False)]
    ut = tab[(True, True)] + tab[(False, True)]
    print()
    print(f"  down preserves => up preserves : {tab[(True, True)]}/{dt}   (failures {tab[(True, False)]})")
    print(f"  up preserves => down preserves : {tab[(True, True)]}/{ut}   (failures {tab[(False, True)]})")
    for k, exs in examples.items():
        print(f"\n  counterexamples for {k}:")
        for e in exs:
            print("   ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

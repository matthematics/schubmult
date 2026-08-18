"""Verify the fiber form of "Gamma is a right ideal", and the derivation of
Lemma forestlift from it.

Setup (writing/schubert_bialgebra_v2.tex).  Taking A = the empty bounded RC
graph of height 1, Definition ringproduct together with Lemma clipempty gives

    R  \sqcupplus  \emptyset_1  =  \sum_{Z(B) = R} B,

the sum over ALL bounded RC graphs B of height h+1 with empty last row and
Z(B) = R (the permutation of B is not fixed: it ranges over the w with
w \downvar{h+1} w_R of equal length).

CLAIM (right ideal, fiber form).  If R_1 =_forest R_2 then for every forest
class Phi,
        #( Z^{-1}(R_1) \cap Phi ) = #( Z^{-1}(R_2) \cap Phi ).

This is exactly the assertion that R_1 \sqcupplus \emptyset_1 - R_2 \sqcupplus
\emptyset_1 lies in Gamma, since Gamma is spanned by differences of forest
equivalent graphs.

DERIVATION of forestlift.  Z is injective on RC_h(w)^0 (Theorem
zero_bijection) and a forest class has a well defined permutation, so each
fiber meets each class at most once; the claim then upgrades to forestlift.
We check the intermediate "at most once" statement too.
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    """height h, permutation perm, empty last row"""
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def rc_all(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h]


def main(N=5):
    # ---- build, for each height h, the full fiber map Z : (height h, empty last row) -> height h-1
    fibers = defaultdict(set)           # R (height h-1)  ->  set of B (height h)
    class_of = {}                       # graph -> forest invariant
    perm_of_class = {}                  # forest invariant -> permutation
    class_wellformed = 0
    class_badperm = 0

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            d = maxd(w)
            if d == 0:
                continue
            for h in range(max(d, 1), n + 1):
                if (w, h) in seen:
                    continue
                seen.add((w, h))
                for B in rc0(w, h):
                    Z = B.zero_out_last_row()
                    fibers[(h, Z)].add(B)
                    phi = B.forest_invariant
                    class_of[B] = phi
                    if phi in perm_of_class:
                        class_wellformed += 1
                        if perm_of_class[phi] != B.perm:
                            class_badperm += 1
                    else:
                        perm_of_class[phi] = B.perm

    # ---- (0) does a forest class determine the permutation?
    print(f"forest class determines permutation: {class_wellformed - class_badperm}/{class_wellformed} (bad {class_badperm})")

    # ---- (1) each fiber meets each forest class at most once
    atmostonce_tot = 0
    atmostonce_bad = 0
    for key, blk in fibers.items():
        c = Counter(class_of[B] for B in blk)
        for phi, m in c.items():
            atmostonce_tot += 1
            if m > 1:
                atmostonce_bad += 1
    print(f"fiber meets each class at most once: {atmostonce_tot - atmostonce_bad}/{atmostonce_tot} (bad {atmostonce_bad})")

    # ---- (2) the right ideal claim: forest equivalent R_1, R_2 have fibers with
    #          the same forest class multiset
    byclass = defaultdict(list)          # (h, forest class of R) -> list of R
    for (h, R) in fibers:
        byclass[(h, R.forest_invariant)].append(R)

    claim_tot = 0
    claim_bad = 0
    examples = []
    for (h, psi), Rs in byclass.items():
        if len(Rs) < 2:
            continue
        sigs = []
        for R in Rs:
            sigs.append(Counter(class_of[B] for B in fibers[(h, R)]))
        claim_tot += 1
        if any(s != sigs[0] for s in sigs):
            claim_bad += 1
            if len(examples) < 4:
                examples.append((h, [tuple(R.perm_word) for R in Rs]))
    print(f"fibers over forest equivalent R have equal class multisets: {claim_tot - claim_bad}/{claim_tot} (bad {claim_bad})")
    for e in examples:
        print("   COUNTEREXAMPLE:", e)

    # ---- (3) the conclusion, checked directly: forestlift
    lift_tot = 0
    lift_bad = 0
    for (h, psi), Rs in byclass.items():
        pool = []
        for R in Rs:
            pool.extend(fibers[(h, R)])
        byperm = defaultdict(list)
        for B in pool:
            byperm[B.perm].append(B)
        for w, blk in byperm.items():
            lift_tot += 1
            if len({class_of[B] for B in blk}) > 1:
                lift_bad += 1
    print(f"forestlift (same perm, Z's forest equivalent => forest equivalent): {lift_tot - lift_bad}/{lift_tot} (bad {lift_bad})")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

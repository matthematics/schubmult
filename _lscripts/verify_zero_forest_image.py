"""Decisive test: what does Z do to a FULL forest class?

For each permutation w and height h with maxd(w) <= h:
  T = RC_h(w)^0  (h rows, empty last row)
  partition T into forest classes  T_Phi  (by P-symbol of the word)
  for each Phi:
     image = Z(T_Phi)
     - how many downstairs forest classes does it meet?
     - is each one met in FULL (i.e. is the image a union of full classes)?

Also recompute the weight generating functions on both sides as a cross-check.
"""

import sys
from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def rc_all(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h]


def main(N=5):
    n_classes = 0
    hist_met = defaultdict(int)      # how many downstairs classes an image meets
    partial = 0                      # image meets a class but not in full
    examples_partial = []
    examples_multi = []
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            d = maxd(w)
            if d == 0:
                continue
            for h in range(d, n + 1):
                T = rc0(w, h)
                if not T:
                    continue
                byclass = defaultdict(list)
                for B in T:
                    byclass[B.forest_invariant].append(B)

                # downstairs universe: all RC graphs of height h-1 over the
                # permutations reachable from w, grouped by forest invariant
                downperms = {B.zero_out_last_row().perm for B in T}
                downuniv = defaultdict(set)
                for wp in downperms:
                    for S in rc_all(wp, h - 1):
                        downuniv[S.forest_invariant].add(S)

                for phi, blk in byclass.items():
                    n_classes += 1
                    image = {B.zero_out_last_row() for B in blk}
                    met = defaultdict(set)
                    for S in image:
                        met[S.forest_invariant].add(S)
                    hist_met[len(met)] += 1
                    if len(met) > 1 and len(examples_multi) < 4:
                        examples_multi.append((tuple(w[:n]), h, len(met), sorted(str(B.perm_word) for B in blk)))
                    for psi, got in met.items():
                        full = downuniv[psi]
                        if got != full:
                            partial += 1
                            if len(examples_partial) < 4:
                                examples_partial.append(
                                    (
                                        tuple(w[:n]),
                                        h,
                                        sorted(str(B.perm_word) for B in blk),
                                        sorted(str(S.perm_word) for S in got),
                                        sorted(str(S.perm_word) for S in full),
                                    )
                                )
    print(f"N={N}: forest classes tested = {n_classes}")
    print(f"       # downstairs classes met, histogram = {dict(sorted(hist_met.items()))}")
    print(f"       image meets a class only partially: {partial}")
    for e in examples_multi:
        print("  MULTI:", e)
    for e in examples_partial:
        print("  PARTIAL:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

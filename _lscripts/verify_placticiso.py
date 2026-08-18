"""Check the conclusion of theorem:placticiso on full Grassmannian RC graphs.

The theorem asserts that on grass_n the squash product is associative and that
P_n : Plac_n -> grass_n is a ring isomorphism.  The step flagged in the proof is

    clp^n(tau_n(T1 (x) T2)) = P_n(T1) box P_n(T2) = P_n(T1 T2),

justified by the rigidity of crystal graph isomorphisms.  We check the two
consequences directly:

  (1) box is associative on full Grassmannian RC graphs of height n;
  (2) box corresponds to the plactic product, i.e. the P-tableau of a squash
      product is the plactic product of the P-tableaux.

A full Grassmannian RC graph of height n is one whose permutation has at most one
descent, at position n.
"""

import sys
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def descents(p):
    return [i for i in range(1, len(p)) if p[i - 1] > p[i]]


def full_grassmannian(n, N):
    out = []
    for k in range(2, N + 1):
        for pw in permutations(range(1, k + 1)):
            p = Permutation(list(pw))
            d = descents(p)
            if len(d) > 1 or (len(d) == 1 and d[0] != n):
                continue
            for R in RCGraph.all_rc_graphs(p, n):
                if len(R) == n:
                    out.append(R)
    # dedupe
    seen, uniq = set(), []
    for R in out:
        if R not in seen:
            seen.add(R)
            uniq.append(R)
    return uniq


def main(n=2, N=4, cap=14):
    G = full_grassmannian(n, N)[:cap]
    print(f"n={n}: {len(G)} full Grassmannian RC graphs of height {n} (capped at {cap})")

    assoc_tot = assoc_ok = 0
    closed = 0
    plac_tot = plac_ok = 0
    bad = []

    for A in G:
        for B in G:
            try:
                AB = A.squash_product(B)
            except Exception:
                continue
            closed += 1
            # (2) plactic compatibility
            try:
                from schubmult import Plactic

                pa = A.p_tableau
                pb = B.p_tableau
                pab = AB.p_tableau
                plac_tot += 1
                prod = Plactic().rs_insert(*(list(pa.row_word) + list(pb.row_word)))
                if tuple(prod.row_word) == tuple(pab.row_word):
                    plac_ok += 1
                elif len(bad) < 3:
                    bad.append(("plactic", tuple(pa.row_word), tuple(pb.row_word),
                                tuple(pab.row_word), tuple(prod.row_word)))
            except Exception as e:
                if len(bad) < 3:
                    bad.append(("plactic-exc", repr(e)))
            for C in G:
                try:
                    left = AB.squash_product(C)
                    right = A.squash_product(B.squash_product(C))
                except Exception:
                    continue
                assoc_tot += 1
                if left == right:
                    assoc_ok += 1
                elif len(bad) < 6:
                    bad.append(("assoc", [tuple(r) for r in A], [tuple(r) for r in B], [tuple(r) for r in C]))

    print(f"  squash product defined on pairs: {closed}")
    print(f"  (1) associativity: {assoc_ok}/{assoc_tot}")
    print(f"  (2) P-tableau of product == plactic product of P-tableaux: {plac_ok}/{plac_tot}")
    for b in bad:
        print("   ", b)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    main(n, N)

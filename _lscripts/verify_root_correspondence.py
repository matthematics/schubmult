"""Does Z induce a root correspondence?

Z preserves the weight vector (Lemma zeroweight), so B and R = Z(B) have the
same number of crossings in every row.  That gives a canonical bijection

    phi:  k-th crossing of row i in B   <->   k-th crossing of row i in R

(within a row, crossings ordered by decreasing column, i.e. word order).

User's suggestion: tracing Algorithm 3 carefully, roots can be matched to the
roots they become, "sort of".  We test:

  (1) does phi preserve the root exactly?  if not, how often, and is the
      failure confined to the crossings touched by the algorithm (the chain
      located in step 7, and the rows inserted into)?
  (2) if phi does not preserve roots, is there SOME bijection that does?
      (compare the multisets of roots of B and of R)
  (3) the application: for an antidiagonal pair R1 -> R2 moving root alpha from
      row i to row i', is the upstairs move phi^{-1}-compatible, i.e. does it
      move the crossing phi^{-1}(alpha-crossing) and land where phi^{-1} says?
"""

import sys
from collections import Counter
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def crossings_by_row(R):
    """row -> list of (i,j) ordered by decreasing column (word order)"""
    out = {}
    for i, row in enumerate(R, start=1):
        out[i] = [(i, L - i + 1) for L in row]
    return out


def crossings(R):
    s = set()
    for i, row in enumerate(R, start=1):
        for L in row:
            s.add((i, L - i + 1))
    return s


def phi_map(B, R):
    """canonical row-order bijection crossings(B) -> crossings(R)"""
    cb, cr = crossings_by_row(B), crossings_by_row(R)
    m = {}
    for i in cr:
        if i not in cb:
            continue
        if len(cb[i]) != len(cr[i]):
            return None
        for a, b in zip(cb[i], cr[i]):
            m[a] = b
    # every crossing of B must be accounted for
    if len(m) != len(crossings(B)):
        return None
    return m


def single_move(c1, c2):
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    return (i, j, ip, jp)


def main(N=5):
    graphs_tested = 0
    phi_defined = 0
    phi_root_exact = 0
    per_crossing_tot = 0
    per_crossing_ok = 0
    multiset_equal = 0
    examples = []

    # application
    app_pairs = 0
    app_phi_ok = 0
    app_src_ok = 0
    app_full_ok = 0

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2 or (w, h) in seen:
                continue
            seen.add((w, h))

            graphs = rc0(w, h)
            if not graphs:
                continue

            Zinv = {}
            for B in graphs:
                R = B.zero_out_last_row()
                Zinv[R] = B
                graphs_tested += 1

                m = phi_map(B, R)
                if m is None:
                    continue
                phi_defined += 1
                allok = True
                for c, d in m.items():
                    per_crossing_tot += 1
                    if B.right_root_at(*c) == R.right_root_at(*d):
                        per_crossing_ok += 1
                    else:
                        allok = False
                if allok:
                    phi_root_exact += 1
                elif len(examples) < 4:
                    bad = [(c, B.right_root_at(*c), d, R.right_root_at(*d))
                           for c, d in sorted(m.items()) if B.right_root_at(*c) != R.right_root_at(*d)]
                    examples.append((tuple(w[:n]), h, [tuple(r) for r in B], [tuple(r) for r in R], bad))

                rb = Counter(B.right_root_at(*c) for c in crossings(B))
                rr = Counter(R.right_root_at(*c) for c in crossings(R))
                if rb == rr:
                    multiset_equal += 1

            # application to antidiagonal pairs
            Rs = list(Zinv)
            for R1 in Rs:
                c1 = crossings(R1)
                for R2 in Rs:
                    if R1 is R2 or R1.perm != R2.perm:
                        continue
                    mv = single_move(c1, crossings(R2))
                    if mv is None:
                        continue
                    i, j, ip, jp = mv
                    if i + j != ip + jp:
                        continue
                    app_pairs += 1
                    B1, B2 = Zinv[R1], Zinv[R2]
                    m1 = phi_map(B1, R1)
                    m2 = phi_map(B2, R2)
                    if m1 is None or m2 is None:
                        continue
                    app_phi_ok += 1
                    inv1 = {v: k for k, v in m1.items()}
                    inv2 = {v: k for k, v in m2.items()}
                    src = inv1.get((i, j))
                    tgt = inv2.get((ip, jp))
                    if src is not None and src[0] == i:
                        app_src_ok += 1
                    if src is not None and tgt is not None:
                        d1, d2 = crossings(B1), crossings(B2)
                        if (d1 - {src}) | {tgt} == d2:
                            app_full_ok += 1

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  graphs tested: {graphs_tested}, row-order bijection defined: {phi_defined}")
    print()
    print(f"  (1) phi preserves ALL roots: {phi_root_exact}/{phi_defined} graphs")
    print(f"      per crossing:            {per_crossing_ok}/{per_crossing_tot}")
    print(f"  (2) root MULTISET of B equals that of R: {multiset_equal}/{phi_defined}")
    print()
    print(f"  (3) antidiagonal pairs: {app_pairs}  (phi defined both sides: {app_phi_ok})")
    print(f"      phi^-1 of the source crossing lies in row i: {app_src_ok}/{app_phi_ok}")
    print(f"      B_2 == B_1 with phi^-1(src) replaced by phi^-1(tgt): {app_full_ok}/{app_phi_ok}")
    for e in examples:
        print(f"\n   ROOT MISMATCH  w={e[0]} h={e[1]}")
        print(f"     B={e[2]}  R={e[3]}")
        for c, rb, d, rr in e[4]:
            print(f"       {c} root {rb}   ->   {d} root {rr}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

"""The Little bump perspective on lemma:forestsimple.

By Proposition little_bump and Proposition zerocrystal, Z is the iteration of
Little bumps, and a Little bump acts LITERALLY ON LETTERS: it decrements the
letter at the bump position (or increments all the others), then, if the result
is not reduced, repeats at the unique position whose deletion restores
reducedness.  Crucially it NEVER REORDERS: letter values change, positions do
not.

USER'S OBSERVATION.  Bumping does not preserve commutation.  But when it fails
to, the PERMUTATION CHANGES: the two outputs still sit at the same two
positions, but they are reduced words for different permutations.

So for U, V reduced words for the same permutation with V = U transposed at
positions p, p+1, and U', V' their bumps, the claim is

    positions:  U' and V' agree off {p, p+1}, and carry the same two letters there
    dichotomy:  perm(U') == perm(V')   <=>   U'_p and U'_{p+1} commute

The contrapositive is exactly what lemma:forestsimple needs: if the two DOWN
words are reduced words for the SAME permutation, then their letters at p, p+1
must commute, so one is the other transposed at p, p+1.

Tested for a single bump and for the full zeroing.
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph

sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from verify_forest_littlebump import bump, bump_position, is_reduced, maxd, perm_of, zero_word


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def transpose_at(w, p):
    return w[:p] + (w[p + 1], w[p]) + w[p + 2 :]


def main(N=5):
    pos_tot = pos_ok = 0
    letters_ok = 0
    dich = Counter()
    dich_bad = []
    z_pos_tot = z_pos_ok = 0
    z_dich = Counter()
    bumppos = Counter()
    ex = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2 or (w, h) in seen:
                continue
            seen.add((w, h))

            graphs = rc0(w, h)
            words = sorted({tuple(B.perm_word) for B in graphs})
            for U in words:
                for p in range(len(U) - 1):
                    if U[p] == U[p + 1]:
                        continue
                    V = transpose_at(U, p)
                    if not is_reduced(V):
                        continue
                    if perm_of(V) != perm_of(U):
                        continue
                    if V not in words:
                        continue
                    if V <= U:
                        continue          # unordered pair, count once

                    # where does the bump start for each?
                    try:
                        bpU, bpV = bump_position(U), bump_position(V)
                        Up, Vp = bump(U), bump(V)
                    except Exception:
                        continue
                    bumppos[(bpU == p or bpU == p + 1, bpU == bpV)] += 1

                    # positions
                    pos_tot += 1
                    off = all(Up[k] == Vp[k] for k in range(len(Up)) if k not in (p, p + 1))
                    if off:
                        pos_ok += 1
                    same_letters = sorted((Up[p], Up[p + 1])) == sorted((Vp[p], Vp[p + 1]))
                    if same_letters:
                        letters_ok += 1

                    # dichotomy
                    commute = abs(Up[p] - Up[p + 1]) >= 2
                    sameperm = perm_of(Up) == perm_of(Vp)
                    dich[(commute, sameperm)] += 1
                    if commute != sameperm and len(dich_bad) < 5:
                        dich_bad.append((tuple(w[:n]), h, U, V, Up, Vp, commute, sameperm))

                    # full zeroing
                    try:
                        Uz, _ = zero_word(U, h)
                        Vz, _ = zero_word(V, h)
                    except Exception:
                        continue
                    z_pos_tot += 1
                    zoff = all(Uz[k] == Vz[k] for k in range(len(Uz)) if k not in (p, p + 1))
                    if zoff:
                        z_pos_ok += 1
                    zcommute = abs(Uz[p] - Uz[p + 1]) >= 2
                    zsame = perm_of(Uz) == perm_of(Vz)
                    z_dich[(zcommute, zsame)] += 1
                    if zcommute != zsame and len(ex) < 5:
                        ex.append((tuple(w[:n]), h, U, V, Uz, Vz, zcommute, zsame))

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  commuting transpositions of a reduced word for the same permutation: {pos_tot}")
    print()
    print("  SINGLE BUMP")
    print(f"    outputs agree off positions p, p+1: {pos_ok}/{pos_tot}")
    print(f"    and carry the same two letters there: {letters_ok}/{pos_tot}")
    print(f"    bump start (at p or p+1?, same for both?): {dict(bumppos)}")
    print(f"    DICHOTOMY (letters commute after bump, permutations agree) -> count:")
    for k in sorted(dich, key=str):
        print(f"       {k}: {dich[k]}")
    print()
    print("  FULL ZEROING")
    print(f"    outputs agree off positions p, p+1: {z_pos_ok}/{z_pos_tot}")
    print(f"    DICHOTOMY -> count:")
    for k in sorted(z_dich, key=str):
        print(f"       {k}: {z_dich[k]}")
    for e in dich_bad:
        print("\n   BUMP DICHOTOMY FAILS:", e)
    for e in ex:
        print("\n   ZERO DICHOTOMY FAILS:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

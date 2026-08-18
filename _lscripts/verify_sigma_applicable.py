"""Does sigma(B_1) even make sense, and is the upstairs move the SAME move?

BUG being corrected: earlier scripts enumerated pairs R1, R2 in the image of Z
WITHOUT requiring w_{R1} = w_{R2}, which Definition antidiag demands.  Here the
pair set is filtered by permutation, and we test the user's objection directly:

  Z changes letters, so even though B_1 has a crossing in row i and B_2 has one
  in row i', they may sit in different columns, i.e. carry a DIFFERENT letter
  than the downstairs move.  In that case the upstairs move is some sigma',
  not sigma.

For each downstairs antidiagonal move sigma = (L, i -> i') we check, upstairs:
  (a) is sigma literally applicable to B_1?  i.e. is (i, L-i+1) in B_1 and is
      (i', L-i'+1) free in B_1?
  (b) is B_2 equal to sigma(B_1) as a set of crossings?
  (c) if B_2 differs from B_1 by a single antidiagonal move sigma', what is it?
      We record (L' - L, i' vs i) to see whether the letter is preserved.
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def crossings(R):
    out = set()
    for i, row in enumerate(R, start=1):
        for L in row:
            out.add((i, L - i + 1))
    return out


def slide_sig(c1, c2):
    """(letter, i, i') if the crossing sets differ by one antidiagonal move."""
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    if i + j != ip + jp:
        return None
    return (i + j - 1, i, ip)


def main(N=5):
    same_perm = 0
    diff_perm = 0
    applicable = 0
    not_applicable = 0
    equals_sigma = 0
    up_sig_kind = Counter()
    letter_shift = Counter()
    examples_bad = []
    examples_diffperm = []
    nontrivial = 0
    nontriv_size = Counter()

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

                Zinv = {}
                for B in rc0(u, m + 1):
                    Zinv[B.zero_out_last_row()] = B
                Rs = list(Zinv)

                for R1 in Rs:
                    c1 = crossings(R1)
                    for R2 in Rs:
                        if R1 is R2:
                            continue
                        sig = slide_sig(c1, crossings(R2))
                        if sig is None:
                            continue
                        L, i, ip = sig
                        if R1.perm != R2.perm:
                            diff_perm += 1
                            if len(examples_diffperm) < 3:
                                examples_diffperm.append(
                                    (tuple(u[:n]), m, sig, tuple(R1.perm[:n]), tuple(R2.perm[:n]))
                                )
                            continue
                        same_perm += 1

                        B1, B2 = Zinv[R1], Zinv[R2]
                        d1, d2 = crossings(B1), crossings(B2)
                        j, jp = L - i + 1, L - ip + 1

                        # how nontrivial is Z here?  compare B1 with R1 directly
                        if d1 != crossings(R1):
                            nontrivial += 1
                            moved = len(d1 ^ crossings(R1))
                            nontriv_size[moved] += 1
                        # does Z change the letter set / the letter L's row?
                        if (i, j) not in crossings(R1):
                            raise RuntimeError("downstairs move not applicable to R1")

                        ok_apply = ((i, j) in d1) and ((ip, jp) not in d1)
                        if ok_apply:
                            applicable += 1
                        else:
                            not_applicable += 1

                        if ok_apply and (d1 - {(i, j)}) | {(ip, jp)} == d2:
                            equals_sigma += 1

                        su = slide_sig(d1, d2)
                        if su is None:
                            up_sig_kind["not a single antidiagonal move"] += 1
                            if len(examples_bad) < 4:
                                examples_bad.append(
                                    (tuple(u[:n]), m, f"down {sig}", "up: multi-crossing",
                                     [tuple(r) for r in B1], [tuple(r) for r in B2])
                                )
                        else:
                            Lu, iu, ipu = su
                            letter_shift[Lu - L] += 1
                            if su == sig:
                                up_sig_kind["identical move"] += 1
                            elif (iu, ipu) == (i, ip):
                                up_sig_kind["same rows, different letter"] += 1
                                if len(examples_bad) < 4:
                                    examples_bad.append(
                                        (tuple(u[:n]), m, f"down {sig}", f"up {su}",
                                         [tuple(r) for r in B1], [tuple(r) for r in B2])
                                    )
                            else:
                                up_sig_kind["different rows"] += 1
                                if len(examples_bad) < 4:
                                    examples_bad.append(
                                        (tuple(u[:n]), m, f"down {sig}", f"up {su}",
                                         [tuple(r) for r in B1], [tuple(r) for r in B2])
                                    )

    print(f"N={N}")
    print(f"  antidiagonal pairs with SAME permutation (the lemma's hypothesis): {same_perm}")
    print(f"  antidiagonal pairs with DIFFERENT permutations (wrongly counted before): {diff_perm}")
    print()
    print(f"  (a) sigma applicable to B_1: {applicable}/{same_perm}   (not applicable {not_applicable})")
    print(f"  (b) B_2 == sigma(B_1) exactly: {equals_sigma}/{same_perm}")
    print(f"  (c) upstairs move kind: {dict(up_sig_kind)}")
    print(f"      letter shift L_up - L_down: {dict(sorted(letter_shift.items()))}")
    print()
    print(f"  nontriviality: Z actually moved crossings in {nontrivial}/{same_perm} instances")
    print(f"      size of symmetric difference B_1 vs R_1: {dict(sorted(nontriv_size.items()))}")
    for e in examples_bad:
        print("   MISMATCH:", e)
    print()
    for e in examples_diffperm:
        print("   different-permutation pair (excluded):", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

"""Classify the FIBER-MATCHED commutation moves only, and test whether the
breakers are chute/ladder moves that lift as themselves.

Framing correction: a commutation neighbor at the level of WORDS does not
determine an RC graph (the compatible sequence is extra data).  The statement
we need only ever compares the canonically matched pair

    B_w in Z^{-1}(R_1),   B'_w in Z^{-1}(R_2)     (same permutation w),

so the classification must be done on those pairs.

A commutation preserves each letter L = i+j-1, so a move that displaces a
single crossing slides it along an ANTI-DIAGONAL: (i,j) -> (i', j') with
i+j = i'+j'.  Displacement +-1 is an elementary chute/ladder move.

We ask, for the fiber-matched pairs:
  (1) how many crossings actually move, downstairs and upstairs?
  (2) crystal vs breaker, by row displacement
  (3) for breakers: is the upstairs move the SAME anti-diagonal displacement
      as the downstairs one?  (the analogue of "B' = same crystal op (B)")
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


def comm_neighbors(word):
    out = []
    for i in range(len(word) - 1):
        if abs(word[i] - word[i + 1]) >= 2:
            out.append(word[:i] + (word[i + 1], word[i]) + word[i + 2 :])
    return out


def crystal_related(R1, R2):
    for i in range(1, len(R1) + 1):
        for op in ("raising_operator", "lowering_operator"):
            try:
                got = getattr(R1, op)(i)
            except Exception:
                got = None
            if got is not None and got == R2:
                return (op, i)
    return None


def move_signature(R1, R2):
    """If exactly one crossing moves, return (letter, i, i'); else None."""
    c1, c2 = crossings(R1), crossings(R2)
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return None
    (i, j), = gone
    (ip, jp), = came
    if i + j != ip + jp:
        return "NONDIAG"
    return (i + j - 1, i, ip)


def main(N=5):
    ncross_down = Counter()
    kind_count = Counter()
    disp = defaultdict(Counter)
    lift_same_disp = Counter()
    lift_ok = Counter()
    ncross_up = defaultdict(Counter)
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
                byword = defaultdict(list)
                for R in Rs:
                    byword[tuple(R.perm_word)].append(R)

                for R1 in Rs:
                    w1 = tuple(R1.perm_word)
                    for w2 in comm_neighbors(w1):
                        for R2 in byword.get(w2, []):
                            if set(fiber[R1]) != set(fiber[R2]):
                                continue
                            sig = move_signature(R1, R2)
                            c1, c2 = crossings(R1), crossings(R2)
                            ncross_down[len(c1 - c2)] += 1
                            if sig is None or sig == "NONDIAG":
                                kind_count["multi-or-nondiag"] += 1
                                continue
                            L, i, ip = sig
                            rel = crystal_related(R1, R2)
                            kind = "crystal" if rel else "breaker"
                            kind_count[kind] += 1
                            disp[kind][ip - i] += 1

                            for w in fiber[R1]:
                                B, Bp = fiber[R1][w], fiber[R2][w]
                                sigu = move_signature(B, Bp)
                                nc = len(crossings(B) - crossings(Bp))
                                ncross_up[kind][nc] += 1
                                lift_ok[kind] += 1
                                if isinstance(sigu, tuple):
                                    Lu, iu, ipu = sigu
                                    if (Lu, iu, ipu) == (L, i, ip):
                                        lift_same_disp[kind] += 1
                                    elif len(examples[kind]) < 4:
                                        examples[kind].append(
                                            (f"down: letter {L} row {i}->{ip}", f"up: letter {Lu} row {iu}->{ipu}")
                                        )
                                elif len(examples[kind]) < 4:
                                    examples[kind].append((f"down: letter {L} row {i}->{ip}", f"up: {sigu}"))

    print("=== fiber-matched pairs: how many crossings move downstairs ===")
    print(f"  {dict(sorted(ncross_down.items()))}")
    print("\n=== classification (single-crossing anti-diagonal moves only) ===")
    print(f"  {dict(kind_count)}")
    print("\n=== row displacement i -> i' ===")
    for k in ("crystal", "breaker"):
        print(f"  {k}: {dict(sorted(disp[k].items()))}")
    print("\n=== upstairs: how many crossings move ===")
    for k in ("crystal", "breaker"):
        print(f"  {k}: {dict(sorted(ncross_up[k].items()))}")
    print("\n=== upstairs move is the IDENTICAL anti-diagonal move ===")
    for k in ("crystal", "breaker"):
        print(f"  {k}: {lift_same_disp[k]}/{lift_ok[k]}")
    print("\n=== samples where upstairs differs from downstairs ===")
    for k, exs in examples.items():
        print(f"  [{k}]")
        for e in exs:
            print("   ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

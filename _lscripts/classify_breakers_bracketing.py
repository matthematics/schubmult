"""Explicit classification of commutation moves into crystal moves and crystal
breakers, in terms of the bracketing of Definition crystal.

For an RCGraph R, row i (1-based) lists LETTERS in decreasing order; the
crossing carrying letter L in row i sits in column L - i + 1.

A single adjacent commutation of word(R) necessarily moves exactly ONE crossing
across a row boundary, preserving its letter (within a row letters strictly
decrease, so a commutation cannot stay inside one row).  We record, for each
such move:

  - the moved letter L, source row i, target row i'
  - the source column j = L-i+1 and target column j' = L-i'+1
  - the bracketing status of the moved crossing in the relevant adjacent pair
    of rows: is its column unpaired?  is it the SMALLEST unpaired one (the one
    f_i is obliged to act on) / the LARGEST unpaired one (for e_i)?
  - whether the move is realized by a crystal operator

Hypothesis to test: the breakers are exactly the moves whose moved crossing is
not the one the crystal operator is obliged to act on (or where the crystal
operator's landing column differs).
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


def bracket(R, i):
    """Bracketing of Definition crystal between rows i and i+1 (1-based).

    Returns (unpaired_in_row_i, unpaired_in_row_i1) as sets of COLUMN indices.
    """
    if i < 1 or i + 1 > len(R):
        return set(), set()
    cols_i = sorted([L - i + 1 for L in R[i - 1]], reverse=True)
    cols_i1 = sorted([L - (i + 1) + 1 for L in R[i]])
    used = set()
    unpaired_i = []
    for b in cols_i:                      # decreasing column
        a = next((a for a in cols_i1 if a > b and a not in used), None)
        if a is None:
            unpaired_i.append(b)
        else:
            used.add(a)
    unpaired_i1 = [a for a in cols_i1 if a not in used]
    return set(unpaired_i), set(unpaired_i1)


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


def main(N=5):
    tab = defaultdict(Counter)
    drow_hist = defaultdict(Counter)
    breaker_reason = Counter()
    examples = defaultdict(list)

    # crystal-move lifting mechanism: is B' literally the crystal image of B?
    mech_tot = mech_ok = 0

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

                graphs = rc_all(u, m)
                byword = defaultdict(list)
                for R in graphs:
                    byword[tuple(R.perm_word)].append(R)

                for R1 in graphs:
                    w1 = tuple(R1.perm_word)
                    c1 = crossings(R1)
                    for w2 in comm_neighbors(w1):
                        for R2 in byword.get(w2, []):
                            c2 = crossings(R2)
                            gone = c1 - c2
                            came = c2 - c1
                            if len(gone) != 1 or len(came) != 1:
                                tab["MALFORMED"][len(gone), len(came)] += 1
                                continue
                            (i, j), = gone
                            (ip, jp), = came
                            L = i + j - 1
                            assert ip + jp - 1 == L, "letter not preserved"

                            rel = crystal_related(R1, R2)
                            kind = "crystal" if rel else "breaker"
                            drow_hist[kind][ip - i] += 1

                            # bracketing status in the adjacent pair the move crosses
                            base = min(i, ip)
                            up_i, up_i1 = bracket(R1, base)
                            if ip > i:      # moving DOWN: f-direction, source row i
                                pool = up_i
                                obliged = min(pool) if pool else None
                                status = (
                                    "unpaired-smallest" if (j in pool and j == obliged)
                                    else "unpaired-not-smallest" if j in pool
                                    else "paired"
                                )
                            else:           # moving UP: e-direction, source row i
                                pool = up_i1
                                obliged = max(pool) if pool else None
                                status = (
                                    "unpaired-largest" if (j in pool and j == obliged)
                                    else "unpaired-not-largest" if j in pool
                                    else "paired"
                                )
                            tab[kind][status] += 1
                            if kind == "breaker":
                                breaker_reason[status, abs(ip - i)] += 1
                                if len(examples[status]) < 3:
                                    examples[status].append(
                                        ([tuple(r) for r in R1], [tuple(r) for r in R2], w1, w2, f"letter {L}: row {i}->{ip}")
                                    )

                            # mechanism check for crystal moves
                            if rel and R1 in fiber and R2 in fiber and set(fiber[R1]) == set(fiber[R2]):
                                op, idx = rel
                                for w in fiber[R1]:
                                    B, Bp = fiber[R1][w], fiber[R2][w]
                                    mech_tot += 1
                                    try:
                                        img = getattr(B, op)(idx)
                                    except Exception:
                                        img = None
                                    if img == Bp:
                                        mech_ok += 1

    print("=== bracketing status of the moved crossing ===")
    for kind in ("crystal", "breaker"):
        print(f"  {kind}: {dict(tab[kind])}")
    print(f"  malformed (symmetric difference not a single crossing): {dict(tab['MALFORMED'])}")
    print("\n=== row displacement i -> i' ===")
    for kind in ("crystal", "breaker"):
        print(f"  {kind}: {dict(sorted(drow_hist[kind].items()))}")
    print("\n=== breaker reasons (status, |row displacement|) ===")
    for k, v in sorted(breaker_reason.items(), key=lambda t: -t[1]):
        print(f"   {k}: {v}")
    print("\n=== crystal-move lifting mechanism: B' == (same crystal op)(B)? ===")
    print(f"   {mech_ok}/{mech_tot}")
    print("\n=== sample breakers by status ===")
    for status, exs in examples.items():
        print(f"  [{status}]")
        for e in exs:
            print("   ", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

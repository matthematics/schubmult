"""Strengthened lockstep test: allow a single braid move as well as a commutation.

Hypothesis: B1, B2 are bounded RC graphs of the same height whose words are
reduced words for the same permutation differing by a single Coxeter move
(a commutation ab <-> ba with |a-b| >= 2, or a braid aba <-> bab with |a-b| = 1),
and whose images under Z have the same permutation.

Claim: the two bumping sequences have equal length and equal intermediate
permutations at every step.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def move_type(u, v):
    """Return 'comm', 'braid', or None according to how u and v differ."""
    if len(u) != len(v) or u == v:
        return None
    diff = [i for i in range(len(u)) if u[i] != v[i]]
    if len(diff) == 2:
        p, q = diff
        if q == p + 1 and u[p] == v[q] and u[q] == v[p] and abs(u[p] - u[q]) >= 2:
            return "comm"
        return None
    if len(diff) == 3:
        p, q, r = diff
        if q == p + 1 and r == p + 2:
            a, b, c = u[p], u[p + 1], u[p + 2]
            d, e, f = v[p], v[p + 1], v[p + 2]
            if a == c and d == f and a == e and b == d and abs(a - b) == 1:
                return "braid"
        return None
    return None


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def bump_sequence(B):
    rc = B
    perms = [rc.perm]
    for _ in range(60):
        if rc.perm.inv == 0 or maxd(rc.perm) != len(rc):
            break
        nxt = rc.little_bump_desc()
        if nxt == rc:
            break
        rc = nxt
        perms.append(rc.perm)
    return perms


def main(n):
    stats = Counter()
    failures = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            reps = {}
            for B in RCGraph.all_rc_graphs(w, h):
                if len(B[-1]) != 0:
                    continue
                reps.setdefault(B.perm_word, B)
            data = {}
            for word, B in reps.items():
                try:
                    data[word] = (B.zero_out_last_row(), bump_sequence(B))
                except Exception:
                    pass
            for u, v in itertools.combinations(data, 2):
                kind = move_type(u, v)
                if kind is None:
                    continue
                img1, seq1 = data[u]
                img2, seq2 = data[v]
                tag = "trivial" if md <= h - 1 else "bump"
                if img1.perm != img2.perm:
                    stats[(tag, kind, "skipped: endpoints differ")] += 1
                    continue
                stats[(tag, kind, "tested")] += 1
                if len(seq1) != len(seq2) or seq1 != seq2:
                    stats[(tag, kind, "FAIL")] += 1
                    if tag == "bump" and len(failures) < 6:
                        failures.append((w, h, kind, u, v, seq1, seq2))
                else:
                    stats[(tag, kind, "ok")] += 1

    print(f"n = {n}")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if failures:
        print("\nsample failures:")
        for w, h, kind, u, v, s1, s2 in failures:
            print(f"  w={w} h={h} move={kind}")
            print(f"    {u} -> {s1}")
            print(f"    {v} -> {s2}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

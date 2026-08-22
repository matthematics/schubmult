"""Key lemma test (permutation form, with the equal-endpoint hypothesis).

Hypothesis: B1, B2 are bounded RC graphs of the same height h whose words lie in
the same commutation class (hence they have the same permutation w), and whose
images under Z have the same permutation.

Claim: the two bumping sequences have the same length and the same intermediate
permutations at every step.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def commutation_class(word):
    """Closure under swaps of adjacent letters differing by at least 2."""
    seen = {tuple(word)}
    stack = [tuple(word)]
    while stack:
        u = stack.pop()
        for i in range(len(u) - 1):
            if abs(u[i] - u[i + 1]) >= 2:
                v = u[: i] + (u[i + 1], u[i]) + u[i + 2 :]
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
    return seen


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def bump_sequence(B):
    """Iterate the single Little bump, recording the permutation at each stage."""
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
            graphs = [B for B in RCGraph.all_rc_graphs(w, h) if len(B[-1]) == 0]
            data = {}
            for B in graphs:
                try:
                    data[B] = (B.zero_out_last_row(), bump_sequence(B))
                except Exception:
                    pass
            classes = {}
            for B in data:
                classes.setdefault(frozenset(commutation_class(B.perm_word)), []).append(B)
            for members in classes.values():
                for B1, B2 in itertools.combinations(members, 2):
                    img1, seq1 = data[B1]
                    img2, seq2 = data[B2]
                    if img1.perm != img2.perm:
                        stats["skipped: endpoints differ"] += 1
                        continue
                    tag = "trivial" if md <= h - 1 else "bump"
                    stats[(tag, "tested")] += 1
                    if len(seq1) != len(seq2):
                        stats[(tag, "FAIL length")] += 1
                        rec = (w, h, B1.perm_word, B2.perm_word, seq1, seq2)
                    elif seq1 != seq2:
                        stats[(tag, "FAIL intermediate perms")] += 1
                        rec = (w, h, B1.perm_word, B2.perm_word, seq1, seq2)
                    else:
                        stats[(tag, "ok")] += 1
                        continue
                    if tag == "bump" and len(failures) < 5:
                        failures.append(rec)

    print(f"n = {n}")
    for k in sorted(stats, key=str):
        print(f"  {k}: {stats[k]}")
    if failures:
        print("\nsample failures:")
        for w, h, w1, w2, s1, s2 in failures:
            print(f"  w={w} h={h}")
            print(f"    word(B1)={w1} -> {s1}")
            print(f"    word(B2)={w2} -> {s2}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

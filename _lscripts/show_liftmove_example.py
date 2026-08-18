"""Print one concrete instance of the lifting statement, for exposition."""

from collections import defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def comm_neighbors(word):
    out = []
    for i in range(len(word) - 1):
        if abs(word[i] - word[i + 1]) >= 2:
            out.append(word[:i] + (word[i + 1], word[i]) + word[i + 2 :])
    return out


def show(R, label):
    print(f"    {label}: rows={[tuple(r) for r in R]}  word={tuple(R.perm_word)}  perm={tuple(R.perm[:6])}")


found = 0
seen = set()
for n in range(2, 6):
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
            byclass = defaultdict(list)
            for R in fiber:
                byclass[R.forest_invariant].append(R)
            for phi, Rs in byclass.items():
                for R1 in Rs:
                    for R2 in Rs:
                        if R1 is R2:
                            continue
                        if tuple(R2.perm_word) not in comm_neighbors(tuple(R1.perm_word)):
                            continue
                        if found >= 2:
                            raise SystemExit
                        found += 1
                        print(f"\n=== instance {found}: u={tuple(u[:n])}, height m={m} ===")
                        print("  DOWNSTAIRS (height m), one commutation apart, same forest class:")
                        show(R1, "R1")
                        show(R2, "R2")
                        print(f"  fiber index set (perms w with w -> u): {[tuple(w[:n+1]) for w in fiber[R1]]}")
                        print("  UPSTAIRS (height m+1), matched by permutation:")
                        for w in fiber[R1]:
                            B, Bp = fiber[R1][w], fiber[R2][w]
                            x, y = tuple(B.perm_word), tuple(Bp.perm_word)
                            pos = [i for i in range(min(len(x), len(y))) if x[i] != y[i]]
                            print(f"    w={tuple(w[:n+1])}")
                            show(B, "  B  (Z(B)=R1) ")
                            show(Bp, "  B' (Z(B')=R2)")
                            print(f"       words differ at positions {pos}: {x} -> {y}")
                            print(f"       same forest class upstairs? {B.forest_invariant == Bp.forest_invariant}")

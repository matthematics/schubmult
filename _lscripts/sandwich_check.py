"""Test the sandwich bounds for N(u,v) against the Tamari-class average.

A = N(u, ddom v)   (bottom of class)      B = N(u, mdom v)   (top of class)
C = A / B <= 1
Ntil = average of N(u,v') over the class.

Claim 1 (user):   (1 - C) * Ntil <= N(u,v) <= (1/C - 1) * Ntil
Claim 2 (interval containment):   C * Ntil <= N(u,v) <= (1/C) * Ntil
"""

import sys
from collections import defaultdict

from schubmult import DSx, Permutation, uncode
from schubmult.symbolic import expand


def mdom(v):
    return ~uncode((~v).theta())


def N(u, v):
    prod = DSx(u) * DSx(v, "z")
    return len({k for k, c in prod.items() if expand(c) != 0})


def main(n):
    perms = list(Permutation.all_permutations(n))

    classes = defaultdict(list)
    for v in perms:
        classes[mdom(v)].append(v)

    # bottom of each class = unique minimal element (lattice congruence class is an interval)
    bots = {}
    for m, cl in classes.items():
        mn = min(cl, key=lambda w: w.inv)
        assert sum(1 for w in cl if w.inv == mn.inv) == 1, f"non-unique min in class {m}"
        bots[m] = mn

    print(f"n={n}: {len(perms)} perms, {len(classes)} Tamari classes")

    bad1lo = bad1hi = bad2lo = bad2hi = 0
    worst1lo = None
    worst1hi = None
    tot = 0
    for u in perms:
        for m, cl in classes.items():
            vals = {v: N(u, v) for v in cl}
            A = vals[bots[m]]
            B = vals[m]
            assert A <= B, f"monotonicity violated u={u} class={m}: A={A} B={B}"
            C = A / B
            Ntil = sum(vals.values()) / len(cl)
            for v, Nv in vals.items():
                tot += 1
                # claim 1
                lo1, hi1 = (1 - C) * Ntil, (1 / C - 1) * Ntil
                if Nv < lo1 - 1e-9:
                    bad1lo += 1
                    r = Nv / lo1 if lo1 > 0 else float("inf")
                    if worst1lo is None or r < worst1lo[0]:
                        worst1lo = (r, u, v, m, bots[m], A, B, Ntil, Nv, len(cl))
                if Nv > hi1 + 1e-9:
                    bad1hi += 1
                    r = Nv / hi1 if hi1 > 0 else float("inf")
                    if worst1hi is None or r > worst1hi[0]:
                        worst1hi = (r, u, v, m, bots[m], A, B, Ntil, Nv, len(cl))
                # claim 2
                if Nv < C * Ntil - 1e-9:
                    bad2lo += 1
                if Nv > Ntil / C + 1e-9:
                    bad2hi += 1

    print(f"  tested {tot} triples (u, class, v)")
    print(f"  claim 1  (1-C)Ntil <= N : {bad1lo} failures")
    print(f"  claim 1  N <= (1/C-1)Ntil: {bad1hi} failures")
    print(f"  claim 2  C*Ntil <= N     : {bad2lo} failures")
    print(f"  claim 2  N <= Ntil/C     : {bad2hi} failures")
    for label, w in (("lower", worst1lo), ("upper", worst1hi)):
        if w:
            r, u, v, m, b, A, B, Ntil, Nv, sz = w
            print(f"  worst claim-1 {label} violation ratio={r:.4f}")
            print(f"     u={u} v={v}  ddom={b} mdom={m} |class|={sz}")
            print(f"     A={A} B={B} C={A / B:.4f} Ntil={Ntil:.4f} N(u,v)={Nv}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)

"""Does x -> (meet(x,a), join(x,a)) inject [e,b] into [e,a] x [a,b] in weak order?

Weak order on S_n is semidistributive; that is what should force injectivity.
Also re-tests |[e,b]| <= |[e,a]| * |[a,b]| for all a <=_L b.
"""

import sys

from schubmult import Permutation


def main(n):
    P = list(Permutation.all_permutations(n))
    leq = {}
    for a in P:
        for b in P:
            leq[(a, b)] = a.weak_order_leq(b)
    meet = {}
    join = {}
    for a in P:
        for b in P:
            meet[(a, b)] = a.weak_order_meet(b)
            join[(a, b)] = a.weak_order_join(b)

    # sanity: meet/join really are the lattice operations
    bad = 0
    for a in P:
        for b in P:
            m, j = meet[(a, b)], join[(a, b)]
            if not (leq[(m, a)] and leq[(m, b)] and leq[(a, j)] and leq[(b, j)]):
                bad += 1
    print(f"S_{n}: {len(P)} perms; meet/join bound-check failures: {bad}")

    # semidistributivity
    sdj = sdm = 0
    for x in P:
        for y in P:
            for z in P:
                if join[(x, y)] == join[(x, z)] and join[(x, y)] != join[(x, meet[(y, z)])]:
                    sdj += 1
                if meet[(x, y)] == meet[(x, z)] and meet[(x, y)] != meet[(x, join[(y, z)])]:
                    sdm += 1
    print(f"  SD-join failures: {sdj}    SD-meet failures: {sdm}")

    # injectivity on [e,b], and the counting consequence
    injbad = 0
    cntbad = 0
    tot = 0
    worst = None
    for a in P:
        for b in P:
            if not leq[(a, b)]:
                continue
            tot += 1
            below_b = [x for x in P if leq[(x, b)]]
            seen = {}
            for x in below_b:
                key = (meet[(x, a)], join[(x, a)])
                if key in seen:
                    injbad += 1
                    if injbad <= 5:
                        print(f"   INJ FAIL a={a} b={b}: {x} and {seen[key]} -> {key}")
                seen[key] = x
            lb = len(below_b)
            la = sum(1 for x in P if leq[(x, a)])
            iv = sum(1 for x in P if leq[(a, x)] and leq[(x, b)])
            if lb > la * iv:
                cntbad += 1
            r = lb / (la * iv)
            if worst is None or r > worst[0]:
                worst = (r, a, b, lb, la, iv)
    print(f"  INJ x->(meet,join) injective on [e,b]: {injbad} failures over {tot} pairs")
    print(f"  |[e,b]| <= |[e,a]|*|[a,b]|: {cntbad} failures")
    r, a, b, lb, la, iv = worst
    print(f"     tightest: a={a} b={b}  {lb} vs {la}*{iv}={la * iv}  ratio={r:.4f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)

"""Tamari-specific: a = ddom v, b = mdom v.
Is x -> (meet(x,a), join(x,a)) injective on [e,b]?
Is Lambda(mdom) <= Lambda(ddom) * m?  Push to n=6.
"""

import sys

from schubmult import Permutation, uncode


def main(n):
    P = list(Permutation.all_permutations(n))
    idx = {w: i for i, w in enumerate(P)}
    inv = {w: frozenset(w.inversion_set) for w in P}

    def leq(a, b):
        return inv[a] <= inv[b]

    def mdom(v):
        return Permutation(list(~uncode((~v).theta())))

    classes = {}
    for v in P:
        classes.setdefault(mdom(v), []).append(v)
    print(f"S_{n}: {len(P)} perms, {len(classes)} classes")

    lam = {}
    for b in set(classes.keys()) | {v for cl in classes.values() for v in cl}:
        lam[b] = sum(1 for x in P if inv[x] <= inv[b])

    bad = 0
    worst = None
    injbad = 0
    for md, cl in classes.items():
        dd = min(cl, key=lambda w: len(inv[w]))
        assert sum(1 for w in cl if len(inv[w]) == len(inv[dd])) == 1
        m = len(cl)
        lhs, rhs = lam[md], lam[dd] * m
        if lhs > rhs:
            bad += 1
            print(f"   FAIL mdom={md} ddom={dd}: {lhs} > {lam[dd]}*{m}")
        r = lhs / rhs
        if worst is None or r > worst[0]:
            worst = (r, md, dd, lhs, lam[dd], m)
        # injectivity of the meet-join map on [e, mdom]
        seen = {}
        for x in P:
            if not leq(x, md):
                continue
            key = (x.weak_order_meet(dd), x.weak_order_join(dd))
            if key in seen:
                injbad += 1
            seen[key] = x

    print(f"  Lambda(mdom) <= Lambda(ddom)*m : {bad} failures / {len(classes)}")
    r, md, dd, lhs, ld, m = worst
    print(f"     tightest: mdom={md} ddom={dd}  {lhs} vs {ld}*{m}={ld * m}  ratio={r:.4f}")
    print(f"  meet-join map injective on [e,mdom]: {injbad} collisions")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

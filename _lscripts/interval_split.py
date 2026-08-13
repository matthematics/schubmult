"""Estimate |[e,mdom v]_L| in terms of |[e,ddom v]_L| and the class size.

Left weak order:  a <=_L b  iff  Inv(a^{-1}) subset Inv(b^{-1}).
Lambda(w) = |[e,w]_L|.

Q1: Lambda(mdom v) <= Lambda(ddom v) * m,  m = |[ddom v, mdom v]_L| ?
Q2: the same for arbitrary a <=_L b:  |[e,b]| <= |[e,a]| * |[a,b]| ?
Q3: is x -> (x /\ a, x \/ a) injective on [e,b]?  (would prove Q2)
"""

import sys
from itertools import permutations

from schubmult import Permutation, uncode


def mdom_of(v, n):
    p = ~uncode((~Permutation(list(v))).theta())
    return tuple(p[i] if i < len(p) else i + 1 for i in range(n))


def inv_mask(w, pairs):
    m = 0
    for b, (i, j) in enumerate(pairs):
        if w[i] > w[j]:
            m |= 1 << b
    return m


def inverse(w):
    n = len(w)
    r = [0] * n
    for i, x in enumerate(w):
        r[x - 1] = i + 1
    return tuple(r)


def main(n):
    P = list(permutations(range(1, n + 1)))
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    # left weak order via inversions of the inverse
    L = {w: inv_mask(inverse(w), pairs) for w in P}

    def leq(a, b):
        return L[a] & ~L[b] == 0

    lam = {}
    for b in P:
        mb = L[b]
        lam[b] = sum(1 for x in P if L[x] & ~mb == 0)

    # Tamari classes
    classes = {}
    for v in P:
        classes.setdefault(mdom_of(v, n), []).append(v)
    bots = {}
    for md, cl in classes.items():
        mn = min(cl, key=lambda w: bin(L[w]).count("1"))
        assert sum(1 for w in cl if bin(L[w]).count("1") == bin(L[mn]).count("1")) == 1
        bots[md] = mn

    print(f"S_{n}: {len(P)} perms, {len(classes)} classes")

    bad = 0
    worst = None
    for md, cl in classes.items():
        dd = bots[md]
        m = len(cl)
        lhs, rhs = lam[md], lam[dd] * m
        if lhs > rhs:
            bad += 1
            if bad <= 5:
                print(f"   FAIL mdom={md} ddom={dd}: {lhs} > {lam[dd]}*{m}={rhs}")
        r = lhs / rhs
        if worst is None or r > worst[0]:
            worst = (r, md, dd, lhs, lam[dd], m)
    r, md, dd, lhs, ld, m = worst
    print(f"  Q1  Lambda(mdom) <= Lambda(ddom)*m : {bad} failures / {len(classes)}")
    print(f"      tightest: mdom={md} ddom={dd}  {lhs} vs {ld}*{m}={ld * m}  ratio={r:.4f}")

    # Q2: arbitrary comparable pairs
    bad2 = 0
    tot2 = 0
    worst2 = None
    for a in P:
        for b in P:
            if not leq(a, b):
                continue
            tot2 += 1
            iv = sum(1 for x in P if leq(a, x) and leq(x, b))
            if lam[b] > lam[a] * iv:
                bad2 += 1
                if bad2 <= 5:
                    print(f"   Q2 FAIL a={a} b={b}: {lam[b]} > {lam[a]}*{iv}")
            r = lam[b] / (lam[a] * iv)
            if worst2 is None or r > worst2[0]:
                worst2 = (r, a, b, lam[b], lam[a], iv)
    r, a, b, lb, la, iv = worst2
    print(f"  Q2  |[e,b]| <= |[e,a]|*|[a,b]| for all a<=_L b : {bad2} failures / {tot2}")
    print(f"      tightest: a={a} b={b}  {lb} vs {la}*{iv}={la * iv}  ratio={r:.4f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

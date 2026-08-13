#!/usr/bin/env python3
"""Is  N_{l(v)}(u, mdom v) <= C * N(u,v) ?   (the paper asserts = O(N(u,v)))

By the double support theorem, the graded piece in degree l(u)+l(v) is
    N_{l(v)}(u,b) = | union_{v'' <=_L b, l(v'')=l(v)}  supp(Sch_u Sch_{v''}) |
    N(u,v)        = | union_{v'' <=_L v}               supp(Sch_u Sch_{v''}) |
Also reported: the raw N(u,b)/N(u,v) ratio and D(v).
"""

import sys

from schubmult import Sx
from schubmult.combinatorics.permutation import Permutation


def D_of(v):
    d, L = 1, 0
    for t in (~v).theta():
        if t:
            d *= t + 1
            L += 1
    return d, L


def covers_down(w):
    wi = ~w
    return [~(wi.swap(d, d + 1)) for d in wi.zero_indexed_descents()]


def below(w):
    seen, fr = {w}, [w]
    while fr:
        nxt = []
        for x in fr:
            for y in covers_down(x):
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        fr = nxt
    return seen


_c = {}


def supp(u, v):
    if (u, v) not in _c:
        _c[(u, v)] = frozenset((Sx(u) * Sx(v)).keys())
    return _c[(u, v)]


def main(n):
    perms = list(Permutation.all_permutations(n))
    bel = {v: below(v) for v in perms}
    worst_graded = (0.0, None)
    worst_full = (0.0, None)
    for v in perms:
        b = v.mul_dominant()
        d, _ = D_of(v)
        lv = v.inv
        graded_src = [x for x in bel[b] if x.inv == lv]
        for u in perms:
            Nuv = len(set().union(*(supp(u, x) for x in bel[v])))
            Ng = len(set().union(*(supp(u, x) for x in graded_src))) if graded_src else 0
            Nb = len(set().union(*(supp(u, x) for x in bel[b])))
            if Ng / Nuv > worst_graded[0]:
                worst_graded = (Ng / Nuv, (u, v, Ng, Nuv, d))
            if Nb / Nuv / d > worst_full[0]:
                worst_full = (Nb / Nuv / d, (u, v, Nb, Nuv, d))
    r, info = worst_graded
    print(f"n={n}:  max N_l(v)(u,b)/N(u,v) = {r:.4f}   at u={info[0]} v={info[1]}  (Ng={info[2]}, N={info[3]}, D={info[4]})")
    r2, i2 = worst_full
    print(f"        max N(u,b)/(D*N(u,v)) = {r2:.4f}   at u={i2[0]} v={i2[1]}  (Nb={i2[2]}, N={i2[3]}, D={i2[4]})")


if __name__ == "__main__":
    for n in range(2, (int(sys.argv[1]) if len(sys.argv) > 1 else 5) + 1):
        main(n)

#!/usr/bin/env python3
"""The commuting-reflections family, where the exponent is extremal.

a_k = s_2 s_4 ... s_{2k} in S_{2k+1}.  These are the class minima realizing the
largest observed log Lambda(b) / log(mean Lambda).  Only the single interval
[e, b]_L is enumerated, so this stays cheap well past the reach of a full S_n
sweep.
"""

import math
import sys

from schubmult.combinatorics.permutation import Permutation


def covers_down(w):
    wi = ~w
    return [~(wi.swap(d, d + 1)) for d in wi.zero_indexed_descents()]


def ideal(b):
    seen = {b}
    frontier = [b]
    while frontier:
        nxt = []
        for w in frontier:
            for x in covers_down(w):
                if x not in seen:
                    seen.add(x)
                    nxt.append(x)
        frontier = nxt
    return seen


def run(k):
    n = 2 * k + 1
    a = Permutation([])
    for j in range(1, k + 1):
        a = a.swap(2 * j - 1, 2 * j)
    b = a.mul_dominant()
    assert a.mul_sortable() == a, "a is not its own class minimum"

    elts = list(ideal(b))
    masks = {w: w.inversion_set for w in elts}
    la = len(ideal(a))
    lb = len(elts)

    amask = masks[a]
    cls = [w for w in elts if amask.issubset(masks[w])]
    for v in cls:
        assert v.mul_sortable() == a and v.mul_dominant() == b, f"class mismatch at {v}"

    lam = {}
    for v in cls:
        mv = masks[v]
        lam[v] = sum(1 for w in elts if masks[w].issubset(mv))

    mean = sum(lam.values()) / len(cls)
    expo = math.log(lb) / math.log(mean)
    print(
        f"k={k} n={n}: Lambda(a)={la:>7} Lambda(b)={lb:>8} |C|={len(cls):>7} "
        f"mean={mean:>10.2f}  Lb/mean={lb / mean:>7.3f}  "
        f"logLb/logmean={expo:.4f}  logLb/logLa={math.log(lb) / math.log(la):.4f}  "
        f"Pr={mean / lb:.4f}",
        flush=True,
    )


if __name__ == "__main__":
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    for k in range(1, hi + 1):
        run(k)

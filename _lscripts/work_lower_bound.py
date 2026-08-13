#!/usr/bin/env python3
"""Empirical lower bound: schubmult's work tracks Lambda(mdom v), not N(u,v).

Family a_k = s_2 s_4 ... s_{2k} in S_{2k+1}:
    N(e, a_k) = Lambda(a_k) = 2^k          (output size)
    Lambda(mdom a_k) = D(a_k) = (2k+1)!!   (dominant hull)
If the measured time tracks (2k+1)!! rather than 2^k, schubmult is not
output-polynomial.
"""

import time

from schubmult import DSx
from schubmult.combinatorics.permutation import uncode


def D_of(v):
    d = 1
    for t in (~v).theta():
        if t:
            d *= t + 1
    return d


print(" k   n    N(e,v)=out    D(v)=Lambda(b)     time(s)    time/out      time/D")
prev = None
for k in range(1, 10):
    n = 2 * k + 1
    v = uncode([0, 1] * k)
    d = D_of(v)
    out_expected = 2**k
    t0 = time.perf_counter()
    prod = DSx([]) * DSx(v, "z")
    t = time.perf_counter() - t0
    got = len(prod)
    print(f" {k}  {n:>3}   {got:>5} ({out_expected})   {d:>10}   {t:>10.4f}  {t / got:>10.2e}  {t / d:>10.2e}")
    if prev:
        print(f"        growth: time x{t / prev[0]:.2f}   out x{got / prev[1]:.2f}   D x{d / prev[2]:.2f}")
    prev = (t, got, d)

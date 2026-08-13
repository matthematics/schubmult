#!/usr/bin/env python3
"""Sampled growth check of  r(n) = max_{u,v} N_{l(v)}(u, mdom v) / N(u,v)."""

import random
import sys

from schubmult import Sx
from schubmult.combinatorics.permutation import Permutation


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


def main(n, samples):
    random.seed(0)
    perms = list(Permutation.all_permutations(n))
    pairs = [(u, v) for u in perms for v in perms]
    if len(pairs) > samples:
        pairs = random.sample(pairs, samples)
    bel = {}
    best = (0.0, None)
    for u, v in pairs:
        for w in (v, v.mul_dominant()):
            if w not in bel:
                bel[w] = below(w)
        b = v.mul_dominant()
        g = [x for x in bel[b] if x.inv == v.inv]
        Nuv = len(set().union(*(supp(u, x) for x in bel[v])))
        Ng = len(set().union(*(supp(u, x) for x in g))) if g else 0
        if Ng / Nuv > best[0]:
            best = (Ng / Nuv, (u, v, Ng, Nuv))
    r, i = best
    print(f"n={n} ({len(pairs)} pairs): max N_l(v)(u,b)/N(u,v) = {r:.4f}   u={i[0]} v={i[1]} Ng={i[2]} N={i[3]}")


if __name__ == "__main__":
    for n in range(4, (int(sys.argv[1]) if len(sys.argv) > 1 else 7) + 1):
        main(n, 4000)

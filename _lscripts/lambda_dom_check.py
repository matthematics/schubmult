#!/usr/bin/env python3
"""Two cheap checks needed for the dominant-hull bound in schubmult.tex.

For sigma in S_n define  Lambda(sigma) = n! / prod_i (sigma(i) - code_i(sigma)).

(1)  sigma dominant  ==>  Lambda(sigma) = Lambda(sigma^{-1}).
(2)  sigma dominant  ==>  Lambda(sigma) <= prod_i (code_i(sigma) + 1).

Pure arithmetic on codes; no interval enumeration.
"""

import sys
from math import factorial, prod


def code_of(arr):
    n = len(arr)
    return [sum(1 for j in range(i + 1, n) if arr[j] < arr[i]) for i in range(n)]


def inverse(arr):
    out = [0] * len(arr)
    for i, a in enumerate(arr):
        out[a - 1] = i + 1
    return out


def big_lambda(arr):
    cd = code_of(arr)
    return factorial(len(arr)) // prod(arr[i] - cd[i] for i in range(len(arr)))


def dominant_from_code(lam, n):
    """The permutation of S_n whose Lehmer code is the partition lam."""
    cd = list(lam) + [0] * (n - len(lam))
    avail = list(range(1, n + 1))
    return [avail.pop(c) for c in cd]


def partitions(n):
    """Partitions that are valid Lehmer codes in S_n, i.e. lam[i] <= n-1-i."""
    out = []

    def rec(pref, maxpart):
        out.append(list(pref))
        i = len(pref)
        if i >= n:
            return
        for p in range(min(maxpart, n - 1 - i), 0, -1):
            rec([*pref, p], p)

    rec([], n - 1)
    return out


def main(hi):
    for n in range(2, hi + 1):
        bad1 = bad2 = tot = 0
        for lam in partitions(n):
            arr = dominant_from_code(lam, n)
            tot += 1
            if big_lambda(arr) != big_lambda(inverse(arr)):
                bad1 += 1
            if big_lambda(arr) > prod(c + 1 for c in lam):
                bad2 += 1
        print(f"n={n:>2}  dominant={tot:>6}   Lambda(s) != Lambda(s^-1): {bad1}   Lambda(s) > prod(theta_i+1): {bad2}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)

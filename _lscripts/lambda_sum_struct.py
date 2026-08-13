"""Probe the structure behind Lambda(alpha+beta) <= Lambda(alpha)Lambda(beta).

P_lambda = {(i,j) : i<j, sigma_lambda(i) < sigma_lambda(j)}  (the forest of non-inversions)
Lambda = number of linear extensions of P_lambda.

Q1: is P_{alpha+beta} = P_alpha  intersect  P_beta  (padded to a common ground set)?
Q2: for "natural" forests (i<j whenever i prec j), is e(P cap Q) <= e(P)e(Q)?
Q3: how does the provable bound prod(a_i+1)prod(b_i+1) compare to Lambda(a)Lambda(b)?
"""

import math
from itertools import permutations

from lambda_sum_bound import Lam, ambient, parts, uncode


def poset(lam, N):
    s = uncode(list(lam) + [0] * (N - len(lam)), N)
    return {(i, j) for i in range(N) for j in range(i + 1, N) if s[i] < s[j]}


def ext_count(rel, N):
    """linear extensions of a relation given as a set of covers/pairs (i<j meaning i before j)"""
    cnt = 0
    for w in permutations(range(N)):
        pos = [0] * N
        for k, x in enumerate(w):
            pos[x] = k
        if all(pos[i] < pos[j] for i, j in rel):
            cnt += 1
    return cnt


def q1(maxpart, maxlen):
    P = parts(maxpart, maxlen)
    bad = 0
    tot = 0
    for a in P:
        for b in P:
            k = max(len(a), len(b))
            aa = list(a) + [0] * (k - len(a))
            bb = list(b) + [0] * (k - len(b))
            s = tuple(x + y for x, y in zip(aa, bb))
            N = max(ambient(s), ambient(a), ambient(b))
            if N > 8:
                continue
            tot += 1
            if poset(s, N) != poset(a, N) & poset(b, N):
                bad += 1
                if bad <= 4:
                    print(f"   a={a} b={b} N={N}")
                    print(f"      P(a+b) = {sorted(poset(s, N))}")
                    print(f"      P(a)&P(b) = {sorted(poset(a, N) & poset(b, N))}")
    print(f"Q1  P(a+b) == P(a) & P(b): {bad} failures out of {tot}")


def q2(N, trials):
    """random natural forests on [N]: parent(i) > i or none, giving i prec parent."""
    import random

    random.seed(0)

    def rand_forest():
        rel = set()
        par = {}
        for i in range(N - 1):
            choices = [None] + list(range(i + 1, N))
            p = random.choice(choices)
            if p is not None:
                par[i] = p
        # transitive closure
        for i in list(par):
            x = i
            seen = set()
            while x in par and x not in seen:
                seen.add(x)
                rel.add((i, par[x]))
                x = par[x]
        return rel

    bad = 0
    for _ in range(trials):
        P, Q = rand_forest(), rand_forest()
        eP, eQ, ePQ = ext_count(P, N), ext_count(Q, N), ext_count(P & Q, N)
        if ePQ > eP * eQ:
            bad += 1
            if bad <= 3:
                print(f"   COUNTEREXAMPLE e(P&Q)={ePQ} > {eP}*{eQ}={eP * eQ}")
                print(f"      P={sorted(P)}  Q={sorted(Q)}")
    print(f"Q2  e(P&Q) <= e(P)e(Q) on random natural forests (N={N}): {bad} failures / {trials}")


def q3(maxpart, maxlen):
    P = parts(maxpart, maxlen)
    worst = None
    for a in P:
        for b in P:
            k = max(len(a), len(b))
            aa = list(a) + [0] * (k - len(a))
            bb = list(b) + [0] * (k - len(b))
            s = tuple(x + y for x, y in zip(aa, bb))
            prodbound = 1
            for x in aa:
                prodbound *= x + 1
            for y in bb:
                prodbound *= y + 1
            L = Lam(s)
            assert L <= prodbound, (a, b, L, prodbound)
            r = prodbound / (Lam(a) * Lam(b))
            if worst is None or r > worst[0]:
                worst = (r, a, b, L, Lam(a) * Lam(b), prodbound)
    r, a, b, L, lp, pb = worst
    print(f"Q3  prod(a_i+1)prod(b_i+1) vs Lambda(a)Lambda(b): worst ratio {r:.3f} at a={a} b={b}")
    print(f"      Lambda(a+b)={L}   Lambda(a)Lambda(b)={lp}   prod bound={pb}")


if __name__ == "__main__":
    q1(3, 3)
    q2(5, 300)
    q3(3, 3)

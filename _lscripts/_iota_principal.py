"""Check: iota(dsch_u^n) = sum_R d_{wt R, u} R. Restricted to R with perm u: is only the principal RC graph present (coef 1)?"""
from itertools import permutations
from schubmult import Permutation, RCGraph, Sx, uncode
from schubmult.symbolic import prod, S
from schubmult.utils.perm_utils import add_perm_dict
import sympy
from collections import Counter

def maxd(w):
    return len(w.trimcode)

g = Sx.genset

def schub_coeff(alpha, u):
    """coefficient of S_u in x^alpha (Schubert expansion)."""
    mono = 1
    for i, a in enumerate(alpha):
        mono *= g[i + 1] ** a
    expansion = Sx(mono)
    return int(expansion.get(u, 0))

stats = Counter(); ex = []
for N in range(2, 6):
    for arr in permutations(range(1, N + 1)):
        u = Permutation(list(arr))
        if u.inv == 0 or u.inv > 5:
            continue
        for n in range(maxd(u), maxd(u) + 2):
            if n > 4:
                continue
            P = RCGraph.principal_rc(u, n)
            for R in RCGraph.all_rc_graphs(u, n):
                alpha = tuple(len(r) for r in R)
                c = schub_coeff(alpha, u)
                stats[("principal", R == P, "coef", c)] += 1
                if R != P and c != 0 and len(ex) < 5:
                    ex.append((tuple(u), n, [tuple(r) for r in R], alpha, c))
print(stats)
for e in ex:
    print("  ", e)

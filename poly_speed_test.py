from time import time
from schubmult import *
from schubmult.symbolic import *
from schubmult.abc import *
from random import shuffle, randint
import sys
from itertools import combinations_with_replacement, permutations
# test time of Schubert polynomial multiplication for random polynomials

perm_sets = {}
def random_multivariate_poly(max_coeff, deg, n_vars):
    arr = list(combinations_with_replacement(list(range(1, n_vars+1)), deg))
    poly_ret = S.Zero
    for term in arr:
        poly_ret += randint(0, max_coeff) * prod([x[i] for i in term])
    return poly_ret

def random_schub_poly(max_coeff, deg, n):
    global perm_sets
    if not perm_sets.get(n):
        perm_sets[n] = {}
        total_perms = [Permutation(perm) for perm in list(permutations(list(range(1, n+1))))]
        for perm in total_perms:
            perm_sets[n][perm.inv] = perm_sets[n].get(perm.inv, []) + [perm]
    perms = perm_sets[n][deg]
    res = S.Zero
    for perm in perms:
        res += randint(0, max_coeff) * Sx(perm)
    return res

n = int(sys.argv[1])
deg = int(sys.argv[2])
max_coeff = int(sys.argv[3])
num = int(sys.argv[4])

a_polys = []
b_polys = []
Sa_polys = []
Sb_polys = []

for i in range(num):
    Sa_polys += [random_schub_poly(max_coeff, deg, n)]
    a_polys += [Sa_polys[-1].as_polynomial().expand()]
    Sb_polys += [random_schub_poly(max_coeff, deg, n)]
    b_polys += [Sb_polys[-1].as_polynomial().expand()]

time_first = time()
for ap, bp in zip(a_polys, b_polys):
    c = expand(ap * bp)
endtime_first = time()

print(f"Regular mult time: {endtime_first - time_first}")

time_first = time()
for Sap, Sbp in zip(Sa_polys, Sb_polys):
    Sa = Sap
    Sb = Sbp
    Sc = (Sa*Sb)
endtime_first = time()

print(f"Schub mult time: {endtime_first - time_first}")


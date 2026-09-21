"""Completeness with a common admissible sigma: for NFs R (height p, maxd<=p) and every sigma in S_p admissible
for R (a_R . v_sigma reduced), is (p, seq, sigma, Q(a_R v_sigma)) injective in R?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from schubmult.utils.perm_utils import is_reduced
from collections import defaultdict

def eg_Q(word):
    P, Q = [], []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return tuple(tuple(r) for r in Q)

def canonical_word(sigma):
    v = []; cur = sigma
    while cur.inv > 0:
        d = max(cur.descents()) + 1
        v.append(d); cur = cur * Permutation.ref_product(d)
    v.reverse()
    return tuple(v)

sig_words = {}
def sigmas(p):
    if p not in sig_words:
        sig_words[p] = [(Permutation(list(s)), canonical_word(Permutation(list(s)))) for s in permutations(range(1, p + 1))]
    return sig_words[p]

tbl = defaultdict(set)
count = 0
for N in range(2, 9):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        p = maxd(w)
        if p < 2 or p > 4:
            continue
        for A in RCGraph.all_rc_graphs(w, p):
            word, seq = A.as_reduced_compatible()
            word = tuple(word)
            for sig, v in sigmas(p):
                full = word + v
                if not is_reduced(list(full)):
                    continue
                count += 1
                tbl[(p, tuple(seq), tuple(sig), eg_Q(full))].add(word)
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("checked", count, "keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:6]:
    print(k, v)

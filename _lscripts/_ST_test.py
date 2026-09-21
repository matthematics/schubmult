"""For each T (rows<=p), S_T = {NF R: seq(R)=seq(T), Inv_ll(R) ⊆ Inv_ll(w_T), a_R v_{sigma_T} reduced,
Q(a_R v_{sigma_T}) = Q(a_T v_{sigma_T})}.  Is |S_T| = 1 (and = {D^p T})?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from schubmult.utils.perm_utils import is_reduced
from collections import defaultdict, Counter

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

def perm_arr(word, n=14):
    w = Permutation.ref_product(*word) if word else Permutation([])
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def inv_ll(word, p):
    a = perm_arr(word)
    return frozenset((i, j) for i in range(p) for j in range(i + 1, p) if a[i] > a[j])

def sigma_of(word, p):
    low = perm_arr(word)[:p]
    order = sorted(range(p), key=lambda i: -low[i])
    return Permutation([order[i] + 1 for i in range(p)])

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

LMAX = 5
# normal forms by (p, seq): enumerate directly all RC graphs of height p with maxd<=p and length<=LMAX
NF = defaultdict(list)
import itertools
def rc_words(p, L):
    # all valid RC graphs with p rows, total L crossings, letters bounded by L+p
    def rows_gen(i, remaining):
        if i > p:
            if remaining == 0:
                yield []
            return
        for k in range(remaining + 1):
            for combo in itertools.combinations(range(i, L + p + 1), k):
                row = tuple(sorted(combo, reverse=True))
                for rest in rows_gen(i + 1, remaining - k):
                    yield [row] + rest
    for rows in rows_gen(1, L):
        word = tuple(x for r in rows for x in r)
        if is_reduced(list(word)):
            w = Permutation.ref_product(*word) if word else Permutation([])
            if maxd(w) <= p:
                seq = tuple(i + 1 for i, r in enumerate(rows) for _ in r)
                NF[(p, seq)].append(word)
for p in range(2, 5):
    for L in range(1, LMAX + 1):
        rc_words(p, L)

stats = Counter(); bad = []
for N in range(2, 9):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > LMAX or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, min(md, 5)):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible()
                word = tuple(word); seq = tuple(seq)
                sig = sigma_of(word, p); v = canonical_word(sig)
                assert is_reduced(list(word + v))
                I = eg_Q(word + v)
                ill = inv_ll(word, p)
                D = tuple(Dp(T, p).as_reduced_compatible()[0])
                S = [R for R in NF[(p, seq)] if inv_ll(R, p) <= ill and is_reduced(list(R + v)) and eg_Q(R + v) == I]
                stats[len(S)] += 1
                if D not in S:
                    print("FAIL", rows[:p], p, tuple(sig), v, D, "Dred", is_reduced(list(D + v)), "Qeq", eg_Q(D + v) == I, "inv", inv_ll(D, p) <= ill, "inNF", D in NF[(p, seq)])
                    raise SystemExit
                if len(S) > 1 and len(bad) < 6:
                    bad.append((rows[:p], p, tuple(sig), D, S))
print(stats)
for b in bad:
    print(b)

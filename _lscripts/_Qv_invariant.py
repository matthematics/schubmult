"""Candidate invariant: I(a) = Q(a . v_w), v_w a canonical reduced word (letters <= p-1) of sigma_w in S_p making the
low values of w.sigma decreasing. Test invariance under p-high bumps and completeness on normal forms."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, reduced_words, maxd
from schubmult.utils.perm_utils import is_reduced
from collections import Counter, defaultdict

def eg_PQ(word):
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
    return tuple(tuple(r) for r in P), tuple(tuple(r) for r in Q)

def perm_arr(word, n):
    w = Permutation.ref_product(*word) if word else Permutation([])
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def canonical_word(sigma):
    """Reduced word (letters as 1-indexed adjacent transpositions) via bubble sort to identity: returns word v
    with product sigma (as functions composed like ref_product)."""
    # find v such that ref_product(*v) == sigma; do it by recursion on right descents (largest first)
    v = []
    cur = sigma
    while cur.inv > 0:
        d = max(cur.descents()) + 1
        v.append(d)
        cur = cur * Permutation.ref_product(d)
    v.reverse()
    # verify
    assert Permutation.ref_product(*v) == sigma
    return tuple(v)

def sigma_w(word, p, n):
    a = perm_arr(word, n)
    low = a[:p]
    order = sorted(range(p), key=lambda i: -low[i])  # indices of low positions in decreasing value order
    # sigma(i) = order[i]+1 ; want (w sigma)(i) = w(sigma(i)) decreasing
    sig = Permutation([order[i] + 1 for i in range(p)])
    return sig

def invariant(word, p, n=12):
    sig = sigma_w(word, p, n)
    v = canonical_word(sig)
    full = tuple(word) + v
    assert is_reduced(list(full)), (word, v)
    P, Q = eg_PQ(full)
    return Q

random.seed(0)
stats = Counter(); ex = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 7:
            continue
        words = list(reduced_words(w))
        if len(words) > 10:
            words = random.sample(words, 10)
        for word in words:
            for p in range(2, maxd(w)):
                I = invariant(word, p)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p:
                        continue
                    nw = ltb_down(word, q)
                    if nw is None:
                        continue
                    I2 = invariant(nw, p)
                    same = I == I2
                    sig_same = sigma_w(word, p, 12) == sigma_w(nw, p, 12)
                    stats[("I invariant", same, "sigma same", sig_same)] += 1
                    if not same and len(ex) < 4:
                        ex.append((word, p, (a, b), nw, I, I2))
print(stats)
for e in ex:
    print(e)

# completeness on normal forms
tbl = defaultdict(set)
for N in range(2, 9):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        p = maxd(w)
        if p < 2 or p > 5:
            continue
        for A in RCGraph.all_rc_graphs(w, p):
            word, seq = A.as_reduced_compatible()
            tbl[(p, tuple(seq), invariant(tuple(word), p))].add(tuple(word))
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("NF completeness: keys", len(tbl), "ambiguous", len(amb))
for k, v in list(amb.items())[:4]:
    print(k, v)

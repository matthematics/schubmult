"""(*) Prefix-locality of the truncated Little map.
L^p(word): iterate word-level Little bumps at the crossing with root (m',m'+1), m'=maxd(current perm), until maxd<=p.
Test: for Y = T u ^p S (rows<=m), prefix_{|a|}(L^p(word(Y))) == L^p(word(T)) == word(D^p T)?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, maxd
from collections import Counter

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])

def Lp(word, p):
    word = tuple(word)
    while True:
        w = perm_of(word)
        m = maxd(w)
        if m <= p:
            return word
        qs = [q for q in range(len(word)) if rroot(word, q) == (m, m + 1)]
        assert len(qs) == 1
        nw = ltb_down(word, qs[0])
        if nw is None:
            return None
        word = nw

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

st = Counter(); ex = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w)
        for Y in RCGraph.all_rc_graphs(w, m):
            rows = [tuple(r) for r in Y]
            for p in range(1, m):
                if all(len(rows[i]) == 0 for i in range(p, m)):
                    continue  # S empty: trivial
                a = tuple(x for r in rows[:p] for x in r)
                if not a:
                    continue
                full = tuple(x for r in rows for x in r)
                L = Lp(full, p)
                if L is None:
                    st["degenerate"] += 1; continue
                T = RCGraph(rows[:p])
                D = tuple(Dp(T, p).as_reduced_compatible()[0])
                La = Lp(a, p)
                st[("L^p(a)==word D^p T", La == D)] += 1
                ok = L[:len(a)] == D
                st[("prefix L^p(as) == D^p T", ok)] += 1
                if not ok and len(ex) < 5:
                    ex.append((rows, p, L, D))
print(st)
for e in ex:
    print("  ", e)

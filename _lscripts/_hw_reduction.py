"""Highest-weight reduction (external Lemma 2). For highest-weight T in RC_{<=p} (e_i T undefined for i<p):
 (a) is D_p(T) determined by (w_T, wt(T))?
 (b) how many highest-weight T per (w, wt, p)?
 (c) is the highest-weight condition preserved by p-high cover bumps?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter, defaultdict

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def is_cover(w, a, b):
    arr = list(w) + list(range(len(w) + 1, b + 2))
    if arr[a - 1] < arr[b - 1]:
        return False
    return not any(arr[b - 1] < arr[c - 1] < arr[a - 1] for c in range(a + 1, b))

tbl = defaultdict(set); cnt = Counter(); pres = Counter(); ex = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                if not is_hw(T, p):
                    continue
                wt = tuple(len(r) for r in rows[:p])
                D = Dp(T, p)
                tbl[(p, tuple(w), wt)].add(tuple(tuple(r) for r in D))
                cnt[(p, tuple(w), wt)] += 1
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p or not is_cover(w, a, b):
                        continue
                    nw = ltb_down(word, q)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    T2 = RCGraph.from_reduced_compatible(list(nw), list(seq))
                    hw2 = is_hw(T2, p)
                    pres[hw2] += 1
                    if not hw2 and len(ex) < 3:
                        ex.append((rows[:p], p, (a, b), [tuple(r) for r in T2]))
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("(a) keys", len(tbl), "ambiguous D_p among hw:", len(amb))
for k, v in list(amb.items())[:4]:
    print("   ", k, v)
print("(b) hw count distribution:", Counter(cnt.values()))
print("(c) hw preserved by p-high cover bumps:", pres)
for e in ex:
    print("   ", e)

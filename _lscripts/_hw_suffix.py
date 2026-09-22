"""For highest-weight T in RC_{<=p}: (i) is the last-descent bump chain a suffix {n-k..n-1} of word positions?
(ii) is word(D^p T) = word(T) - v with v weakly increasing along positions (i.e. a 'staircase' decrement)?
(iii) is v determined by (p, wt, w)?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd, rroot, ltb_down
from collections import Counter, defaultdict

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

st = Counter(); ex = []; vt = defaultdict(set)
for N in range(2, 9):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w)
        for p in range(2, m):
            for T in RCGraph.all_rc_graphs(w, m):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, m)) or not is_hw(T, p):
                    continue
                word = tuple(T.as_reduced_compatible()[0]); n = len(word)
                q = next(k for k in range(n) if rroot(word, k) == (m, m + 1))
                nw, chain = ltb_down(word, q, return_chain=True)
                if nw is not None:
                    suffix = sorted(chain) == list(range(n - len(chain), n))
                    st[("chain is suffix", suffix)] += 1
                    st[("chain starts at last position", q == n - 1)] += 1
                    if not suffix and len(ex) < 4:
                        ex.append((rows[:p], word, q, chain))
                D = Dp(T, p)
                dw = tuple(D.as_reduced_compatible()[0])
                v = tuple(a - b for a, b in zip(word, dw))
                st[("v weakly increasing", all(v[i] <= v[i + 1] for i in range(n - 1)))] += 1
                st[("v nonneg", all(x >= 0 for x in v))] += 1
                vt[(p, tuple(len(r) for r in rows[:p]), tuple(w))].add(v)
amb = sum(1 for s in vt.values() if len(s) > 1)
print(st)
print("v determined by (p,wt,w)? keys", len(vt), "ambiguous", amb)
for e in ex:
    print("  ", e)

"""Lemma C' : for hw T in RC_{<=p}, #{hw A in RC_p : trim A = trim D^p T, wt A = wt T, extwt A = extwt T, w_A in Tcal(w_T)} == 1 ?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations, combinations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import is_reduced
from _lb_common import maxd
from collections import defaultdict, Counter

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def extwt(T, p):
    seen = {T}; stack = [T]; best = tuple(len(r) for r in T)[:p]
    while stack:
        X = stack.pop()
        wt = tuple(len(r) for r in X)[:p]
        if wt < best:
            best = wt
        for i in range(1, p):
            Y = X.lowering_operator(i)
            if Y is not None and Y not in seen:
                seen.add(Y); stack.append(Y)
    return best

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def transition_set(z, n):
    if maxd(z) < n:
        return {z}
    out = set(); stack = [z]; seen = set()
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u)
        if maxd(u) < n:
            out.add(u); continue
        au = arr(u)
        s = max(j for j in range(n + 1, len(au) + 1) if au[j - 1] < au[n - 1])
        b = au[:]; b[n - 1], b[s - 1] = b[s - 1], b[n - 1]
        uts = Permutation(b)
        for i in range(1, n):
            c = b[:]; c[i - 1], c[n - 1] = c[n - 1], c[i - 1]
            up = Permutation(c)
            if up.inv == uts.inv + 1:
                stack.append(up)
    return out

def Tcal(w, p):
    cur = {w}
    for n in range(maxd(w), p, -1):
        nxt = set()
        for u in cur:
            nxt |= transition_set(u, n)
        cur = nxt
    return cur

# index all hw NFs in RC_p by (p, trim rows, wt1)
def rc_graphs(p, L):
    def rows_gen(i, remaining):
        if i > p:
            if remaining == 0:
                yield []
            return
        for k in range(remaining + 1):
            for combo in combinations(range(i, L + p + 1), k):
                row = tuple(sorted(combo, reverse=True))
                for rest in rows_gen(i + 1, remaining - k):
                    yield [row] + rest
    for rows in rows_gen(1, L):
        word = [x for r in rows for x in r]
        if is_reduced(word):
            T = RCGraph(rows)
            if maxd(T.perm) <= p:
                yield T

LMAX = 6
NF = defaultdict(list)
for p in range(2, 5):
    for L in range(1, LMAX + 1):
        for A in rc_graphs(p, L):
            if is_hw(A, p):
                rows = [tuple(r) for r in A]
                NF[(p, tuple(rows[1:]), len(rows[0]))].append((A, extwt(A, p)))

stats = Counter(); ex = []
for N in range(2, 9):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > LMAX or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, min(md, 5)):
            Tset = None
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)) or not is_hw(T, p):
                    continue
                if Tset is None:
                    Tset = Tcal(w, p)
                D = Dp(T, p); drows = [tuple(r) for r in D]
                e = extwt(T, p)
                cands = [A for A, eA in NF[(p, tuple(drows[1:]), len(drows[0]))] if eA == e and A.perm in Tset]
                stats[len(cands)] += 1
                if D not in cands:
                    stats["D missing (NF list incomplete)"] += 1
                if len(cands) > 1 and len(ex) < 5:
                    ex.append((rows[:p], p, arr(w, 8), drows, [[tuple(r) for r in A] for A in cands]))
print(stats)
for x in ex:
    print("  ", x)

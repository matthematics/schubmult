"""Candidate Lemma C and refinements.
For T in RC_{<=p} (height h=maxd>p): 𝒯(w_T) = endpoints of iterated last-descent transitions T_h,...,T_{p+1}.
v = w(trim D_p T), a = wt_1(T). Count w' in 𝒯(w_T) with v in P_1(w') and ell(w')-ell(v)=a.
Also: count over all w' in 𝒯(w_T) (any a) ; and the set 𝒯(w_T) vs 𝒯(w_{T'}) for T'=(Z Y)_{<=p}."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations, combinations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from collections import Counter

def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def transition_set(z, n):
    """T_n(z): permutations with maxd<n reachable by transition steps at n (same length)."""
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
        uts = u * Permutation.ref_product(*[]) if False else None
        # u t_{ns}
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

dec_cache = {}
def decreasing_products(a, M=12):
    if (a, M) not in dec_cache:
        S = set()
        for comb in combinations(range(1, M + 1), a):
            S.add(Permutation.ref_product(*sorted(comb, reverse=True)) if a else Permutation([]))
        dec_cache[(a, M)] = S
    return dec_cache[(a, M)]

def shift_up(v):
    return Permutation([1] + [x + 1 for x in v])

def in_P1(wp, v, a):
    c = wp * (~shift_up(v))
    return c.inv == a and c in decreasing_products(a)

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def trim(X):
    return RCGraph([tuple(x - 1 for x in r) for r in list(X)[1:]])

stats = Counter(); ex = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            Tset = Tcal(w, p)
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                D = Dp(T, p)
                v = trim(D).perm
                a = len(rows[0])
                assert D.perm in Tset
                cands = [wp for wp in Tset if in_P1(wp, v, a)]
                stats[("#cands with v in P1, length a", len(cands))] += 1
                if len(cands) > 1 and len(ex) < 5:
                    ex.append((rows[:p], p, arr(w, 8), arr(v, 6), a, [arr(x, 8) for x in cands], arr(D.perm, 8)))
print(stats)
for e in ex:
    print(e)

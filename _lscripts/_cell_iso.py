"""User's isomorphism: T (perm u), T' (perm v) highest weight with isomorphic crystals (same shape, same extwt);
sigma: root at cell c of T -> root at cell c of T'.  Test 'commutes with Little bumps':
for each cell c whose root r is a cover of u and sigma(r) a cover of v (both p-high or any): apply bumps at those roots;
(1) are both results hw with isomorphic crystals? (2) same set of decremented cells? (3) do the new root tableaux
again correspond cell-to-cell consistently (i.e. cells unchanged by the bump keep sigma-related roots)?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter, defaultdict

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])
def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]
def is_cover(w, a, b):
    ar = arr(w)
    if ar[a - 1] < ar[b - 1]:
        return False
    return not any(ar[b - 1] < ar[c - 1] < ar[a - 1] for c in range(a + 1, b))
def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))
def extwt(T, p):
    seen = {T}; stack = [T]; best = tuple(len(r) for r in T)[:p]
    while stack:
        X = stack.pop(); wt = tuple(len(r) for r in X)[:p]
        if wt < best: best = wt
        for i in range(1, p):
            Y = X.lowering_operator(i)
            if Y is not None and Y not in seen:
                seen.add(Y); stack.append(Y)
    return best, len(seen)

def eg_Q(word):
    P, Q = [], []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]; bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x: x = y
            else: row[k] = x; x = y
            r += 1
    return Q

def cells(word):
    """cell -> position (HY convention: insert reversed word)."""
    Q = eg_Q(tuple(reversed(word))); n = len(word)
    return {(i, j): n - t for i, row in enumerate(Q) for j, t in enumerate(row)}

def root_tab(word):
    return {c: rroot(word, q) for c, q in cells(word).items()}

def bump_cells(word, q):
    nw, ch = ltb_down(word, q, return_chain=True)
    if nw is None: return None, None
    pos2cell = {q_: c for c, q_ in cells(word).items()}
    return nw, frozenset(pos2cell[k] for k in ch)

# collect hw graphs by (p, wt, extwt)
groups = defaultdict(list)
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 6 or w.inv == 0: continue
        md = maxd(w)
        for p in range(2, md + 1):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)) or not is_hw(T, p): continue
                e, sz = extwt(T, p)
                word, seq = T.as_reduced_compatible()
                groups[(p, tuple(len(r) for r in rows[:p]), e, sz)].append((tuple(word), tuple(seq), w))

stats = Counter(); ex = []
for key, lst in groups.items():
    p = key[0]
    for i in range(len(lst)):
        for j in range(len(lst)):
            if i == j: continue
            (wa, sa, u), (wb, sb, v) = lst[i], lst[j]
            RA, RB = root_tab(wa), root_tab(wb)
            if set(RA) != set(RB): stats["different cell sets"] += 1; continue
            cA, cB = cells(wa), cells(wb)
            for c in RA:
                ra, rb = RA[c], RB[c]
                covA = ra[0] > p and is_cover(u, *ra); covB = rb[0] > p and is_cover(v, *rb)
                stats[("cover status agree", covA == covB)] += 1
                if not (covA and covB): continue
                na, SA = bump_cells(wa, cA[c]); nb, SB = bump_cells(wb, cB[c])
                if na is None or nb is None or not valid_rc(na, sa) or not valid_rc(nb, sb):
                    stats["degenerate"] += 1; continue
                TA = RCGraph.from_reduced_compatible(list(na), list(sa)); TB = RCGraph.from_reduced_compatible(list(nb), list(sb))
                hwA, hwB = is_hw(TA, p), is_hw(TB, p)
                stats[("both hw", hwA and hwB)] += 1
                if not (hwA and hwB): continue
                eA, eB = extwt(TA, p), extwt(TB, p)
                stats[("results isomorphic (extwt,size)", eA == eB)] += 1
                stats[("same decremented cells", SA == SB)] += 1
                if SA != SB and len(ex) < 4:
                    ex.append((wa, wb, c, ra, rb, sorted(SA), sorted(SB)))
                # consistency: cells outside SA∪SB keep the same pairing? i.e. new root tableaux define sigma' extending sigma on untouched cells
                RA2, RB2 = root_tab(na), root_tab(nb)
                same_outside = all((RA2[d] == RA[d]) == (RB2[d] == RB[d]) for d in RA)
                stats[("cells changed agree", same_outside)] += 1
for k, v in sorted(stats.items(), key=lambda x: str(x[0])):
    print(v, k)
for e in ex: print("  ", e)

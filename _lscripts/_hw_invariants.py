"""On highest-weight T in RC_{<=p}: test candidate invariants I(T) for
 (inv) invariance under p-high cover bumps (which preserve hw), and
 (det) whether (p, wt, I) determines D^p(T)."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter, defaultdict

def eg_insert(word):
    P = []
    for x in word:
        r = 0
        while True:
            if r == len(P):
                P.append([x]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return tuple(tuple(r) for r in P)

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

def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def std(vals):
    s = sorted(vals); return tuple(s.index(v) + 1 for v in vals)

def invariants(T, p):
    w = T.perm; a = arr(w)
    word = tuple(T.as_reduced_compatible()[0])
    P = eg_insert(tuple(reversed(word)))
    low = a[:p]
    out = {}
    out["std_low"] = std(low)
    out["fixed_low"] = tuple(i for i in range(p) if a[i] == i + 1)
    out["inv_ll"] = frozenset((i, j) for i in range(p) for j in range(i + 1, p) if a[i] > a[j])
    out["low_vals_le_p"] = tuple(v if v <= p else 0 for v in low)
    out["P_le_p_cells"] = frozenset((i, j, v) for i, row in enumerate(P) for j, v in enumerate(row) if v <= p)
    out["P_lt_p_cells"] = frozenset((i, j, v) for i, row in enumerate(P) for j, v in enumerate(row) if v < p)
    out["P_le_p_shape"] = tuple(sum(1 for v in row if v <= p) for row in P)
    out["winv_low_pos"] = tuple(a.index(v) + 1 if v in a else 0 for v in range(1, p + 1))  # positions of values 1..p
    out["winv_low_pos_std"] = std([a.index(v) + 1 for v in range(1, p + 1)])
    return out

inv_stats = Counter(); det = defaultdict(lambda: defaultdict(set))
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)) or not is_hw(T, p):
                    continue
                I = invariants(T, p)
                wt = tuple(len(r) for r in rows[:p])
                D = tuple(tuple(r) for r in Dp(T, p))
                for k, v in I.items():
                    det[k][(p, wt, v)].add(D)
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p or not is_cover(w, a, b):
                        continue
                    nw = ltb_down(word, q)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    T2 = RCGraph.from_reduced_compatible(list(nw), list(seq))
                    I2 = invariants(T2, p)
                    for k in I:
                        inv_stats[(k, I[k] == I2[k])] += 1
print("INVARIANCE under p-high cover bumps (hw):")
for k in sorted({k for k, _ in inv_stats}):
    print(f"  {k}: same={inv_stats[(k, True)]}, changed={inv_stats[(k, False)]}")
print("DETERMINATION of D^p by (p, wt, I):")
for k, tb in det.items():
    amb = sum(1 for v in tb.values() if len(v) > 1)
    print(f"  {k}: keys={len(tb)}, ambiguous={amb}")

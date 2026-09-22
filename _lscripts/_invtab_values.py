"""In InversionsTableau coordinates: under a p-high cover bump, do all changed (root,label) entries share a
common coordinate (the bumped wire)? And is the set of roots not containing that wire, with labels, unchanged?
Then along D^p: which value-coordinate sub-tableau is preserved?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.combinatorics.inversions_tableau import InversionsTableau
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter

def inv_tab(T):
    return dict(InversionsTableau.from_rc_graph(T)._dict)

def is_cover(w, a, b):
    arr = list(w) + list(range(len(w) + 1, b + 2))
    if arr[a - 1] < arr[b - 1]:
        return False
    return not any(arr[b - 1] < arr[c - 1] < arr[a - 1] for c in range(a + 1, b))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

T0 = RCGraph([(2,), (4,)])
print("example coords: perm", T0.perm, "invtab", inv_tab(T0))
st = Counter(); wires = Counter(); st2 = Counter()
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        aw = list(w) + list(range(len(w) + 1, 12))
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                I = inv_tab(T)
                D = Dp(T, p); ID = inv_tab(D)
                # value-coordinate: roots (x,y) values; low values = values at positions<=p? test sub-tableau on roots with both values <= p
                for name, pred in (("both<=p", lambda r: r[0] <= p and r[1] <= p), ("min<=p", lambda r: min(r) <= p), ("max<=p", lambda r: max(r) <= p)):
                    sub = {r: v for r, v in I.items() if pred(r)}; subD = {r: v for r, v in ID.items() if pred(r)}
                    st2[(name, "equal", sub == subD)] += 1
                    st2[(name, "D subset T", all(I.get(r) == v for r, v in subD.items()))] += 1
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p or not is_cover(w, a, b):
                        continue
                    nw = ltb_down(word, q)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    T2 = RCGraph.from_reduced_compatible(list(nw), list(seq))
                    I2 = inv_tab(T2)
                    changed = [r for r in set(I) | set(I2) if I.get(r) != I2.get(r)]
                    common = set(changed[0]) if changed else set()
                    for r in changed[1:]:
                        common &= set(r)
                    st[("changed share a wire", len(common) >= 1)] += 1
                    if len(common) == 1:
                        x = next(iter(common))
                        # relate x to a (position) : x == w(a)? or w(b)?
                        wires[("x==w(a)", x == aw[a - 1], "x==w(b)", x == aw[b - 1], "x==a", x == a)] += 1
print(st)
print(wires)
print(st2)

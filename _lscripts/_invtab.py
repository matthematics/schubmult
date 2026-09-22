"""Inversions tableau (root -> row label) under p-high cover bumps and under D^p.
(1) For a bump T->T2 at (a,b): which entries of the inversions tableau change?
(2) Along D^p: which entries survive from T to D^p(T)? Is the sub-tableau on roots (i,j) with j<=p preserved?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.combinatorics.inversions_tableau import InversionsTableau
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter

def inv_tab(T):
    return dict(InversionsTableau.from_rc_graph(T)._dict)

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

st1 = Counter(); st2 = Counter(); ex = []
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
                I = inv_tab(T)
                # (2) along D^p
                D = Dp(T, p); ID = inv_tab(D)
                ll = {r: v for r, v in I.items() if r[1] <= p}
                llD = {r: v for r, v in ID.items() if r[1] <= p}
                st2[("low-low sub-tableau: D subset of T", all(llD.get(r) == v for r, v in llD.items() if r in ll) and set(llD) <= set(ll))] += 1
                st2[("low-low equal", ll == llD)] += 1
                # roots (i,j) with i<=p: labels multiset per i
                lab_i = tuple(sorted(Counter(v for r, v in I.items() if r[0] == i).items()) for i in range(1, p + 1))
                lab_iD = tuple(sorted(Counter(v for r, v in ID.items() if r[0] == i).items()) for i in range(1, p + 1))
                st2[("labels by first coord<=p multiset", lab_i == lab_iD)] += 1
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
                    changed = {r for r in set(I) | set(I2) if I.get(r) != I2.get(r)}
                    st1[("changed roots all contain a", all(a in r for r in changed))] += 1
                    st1[("low-low unchanged", all(I.get(r) == I2.get(r) for r in set(I) | set(I2) if r[1] <= p))] += 1
                    # removed (a,b), added (c,a); relabeling pattern
                    removed = set(I) - set(I2); added = set(I2) - set(I)
                    st1[("removed==(a,b)", removed == {(a, b)})] += 1
                    st1[("added=(c,a) c<a", len(added) == 1 and next(iter(added))[1] == a)] += 1
                    if not all(a in r for r in changed) and len(ex) < 3:
                        ex.append((rows[:p], p, (a, b), sorted(changed)))
print("(1) bumps:", st1)
for e in ex:
    print("   ", e)
print("(2) along D^p:", st2)

"""Left-root inversions tableau (root transported through the prefix) under bumps and along D^p."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter

def left_tab(T):
    d = {}
    for i in range(T.perm.inv):
        r = T.left_to_right_left_inversion(i)
        d[(min(r), max(r))] = T.left_to_right_inversion_coords(i)[0]
    return d

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
print("example: perm", T0.perm, "left tab", left_tab(T0))
st = Counter(); st2 = Counter(); ex = []
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
                I = left_tab(T)
                D = Dp(T, p); ID = left_tab(D)
                for name, pred in (("both<=p", lambda r: r[1] <= p), ("min<=p", lambda r: r[0] <= p), ("max>p", lambda r: r[1] > p)):
                    sub = {r: v for r, v in I.items() if pred(r)}; subD = {r: v for r, v in ID.items() if pred(r)}
                    st2[(name, "equal", sub == subD)] += 1
                    st2[(name, "D subset T", all(I.get(r) == v for r, v in subD.items()))] += 1
                    st2[(name, "T subset D", all(ID.get(r) == v for r, v in sub.items()))] += 1
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    if a <= p or not is_cover(w, a, b):
                        continue
                    nw = ltb_down(word, q)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    T2 = RCGraph.from_reduced_compatible(list(nw), list(seq))
                    I2 = left_tab(T2)
                    changed = [r for r in set(I) | set(I2) if I.get(r) != I2.get(r)]
                    common = set(changed[0]) if changed else set()
                    for r in changed[1:]:
                        common &= set(r)
                    st[("changed share a wire", len(common) >= 1)] += 1
                    if len(common) == 1:
                        x = next(iter(common))
                        st[("wire", "w(a)" if x == aw[a - 1] else "w(b)" if x == aw[b - 1] else "a" if x == a else "other")] += 1
                        # do the unchanged-wire roots keep labels; do roots on wire x keep the multiset of labels?
                        lx = sorted(v for r, v in I.items() if x in r); lx2 = sorted(v for r, v in I2.items() if x in r)
                        st[("labels on wire x multiset same", lx == lx2)] += 1
                    elif len(ex) < 3:
                        ex.append((rows[:p], p, (a, b), changed))
print(st)
for e in ex:
    print("  ", e)
print(st2)

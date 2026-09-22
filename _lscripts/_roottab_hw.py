"""Root tableau (cell -> right root) of highest-weight T, before/after one last-descent bump (lmap)
and after full D^p. Print small examples to see the cellwise relation."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.combinatorics.root_tableau import RootTableau
from _lb_common import maxd, rroot

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def eg_Q(word):
    P, Q = [], []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return Q

def root_grid(T):
    word = tuple(T.as_reduced_compatible()[0]); n = len(word)
    Q = eg_Q(tuple(reversed(word)))
    g = {}
    for i, row in enumerate(Q):
        for j, t in enumerate(row):
            g[(i, j)] = rroot(word, n - t)
    return g

def show(g):
    rows = {}
    for (i, j), r in g.items():
        rows.setdefault(i, {})[j] = r
    return " | ".join(" ".join(f"{rows[i][j][0]}{rows[i][j][1]}" for j in sorted(rows[i])) for i in sorted(rows))

def arr(w, n=8):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return "".join(map(str, a[:n]))

count = 0
for N in range(3, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 5 or w.inv < 2:
            continue
        md = maxd(w)
        for p in (2, 3):
            if md <= p:
                continue
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)) or not is_hw(T, p):
                    continue
                T1 = T.little_bump_desc()
                D = Dp(T, p)
                print(f"p={p} w={arr(w)} T={rows[:p]}  roots: {show(root_grid(T))}")
                print(f"      lmap: w={arr(T1.perm)} T1={[tuple(r) for r in T1][:p]}  roots: {show(root_grid(T1))}")
                print(f"      D^p : w={arr(D.perm)} D={[tuple(r) for r in D]}  roots: {show(root_grid(D))}")
                count += 1
                if count >= 28:
                    raise SystemExit

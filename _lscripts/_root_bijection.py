"""Test: for hw T (perm u), T' (perm v) in RC_p with same weight and extwt (isomorphic Demazure crystals),
is the crystal isomorphism phi induced by a single bijection sigma: I(u)->I(v) of inversion sets,
i.e. invtab(phi(R)) = invtab(R) o sigma^{-1} for all R in the component?"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.combinatorics.inversions_tableau import InversionsTableau
from _lb_common import maxd
from collections import Counter, defaultdict

def is_hw(T, p):
    return all(T.raising_operator(i) is None for i in range(1, p))

def invtab(T):
    return dict(InversionsTableau.from_rc_graph(T)._dict)

def component(T, p):
    """dict path(tuple of i's) -> element, BFS by lowering operators; also extwt."""
    comp = {(): T}; order = [()]
    seen = {T: ()}
    k = 0
    while k < len(order):
        path = order[k]; X = comp[path]; k += 1
        for i in range(1, p):
            Y = X.lowering_operator(i)
            if Y is not None and Y not in seen:
                seen[Y] = path + (i,); comp[path + (i,)] = Y; order.append(path + (i,))
    ext = min(tuple(len(r) for r in X)[:p] for X in comp.values())
    return comp, ext, seen

# group hw by (p, wt, ext)
groups = defaultdict(list)
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 6 or w.inv == 0:
            continue
        p = maxd(w)
        if p < 2 or p > 4:
            continue
        for T in RCGraph.all_rc_graphs(w, p):
            if not is_hw(T, p):
                continue
            comp, ext, seen = component(T, p)
            groups[(p, tuple(len(r) for r in T), ext)].append((T, comp, seen))

stats = Counter(); ex = []
for key, lst in groups.items():
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            T, comp, seen = lst[i]; T2, comp2, seen2 = lst[j]
            # match elements by lowering-path (crystal iso must respect f_i); need same path sets
            if set(seen.values()) != set(seen2.values()):
                # paths differ (BFS first-path); match via canonical: for each element of comp, follow its path in comp2
                pass
            phi = {}
            ok = True
            for X, path in seen.items():
                Y = T2
                for a in path:
                    Y = Y.lowering_operator(a)
                    if Y is None:
                        ok = False; break
                if not ok:
                    break
                phi[X] = Y
            if not ok or len(set(phi.values())) != len(comp2):
                stats["not isomorphic via paths"] += 1
                continue
            # find sigma: for each root r of u, its label in every R must equal label of sigma(r) in phi(R)
            Iu = sorted(invtab(T)); Iv = sorted(invtab(T2))
            # signature of root r: tuple of labels across component elements in fixed order
            elems = list(seen.keys())
            sig_u = {r: tuple(invtab(X)[r] for X in elems) for r in Iu}
            sig_v = {r: tuple(invtab(phi[X])[r] for X in elems) for r in Iv}
            inv_v = defaultdict(list)
            for r, s in sig_v.items():
                inv_v[s].append(r)
            unique = all(len(inv_v[sig_u[r]]) == 1 for r in Iu if sig_u[r] in inv_v)
            exists = all(sig_u[r] in inv_v for r in Iu)
            stats[("sigma exists", exists, "unique", exists and unique, "same perm", T.perm == T2.perm)] += 1
            if not exists and len(ex) < 3:
                ex.append((key, [tuple(r) for r in T], [tuple(r) for r in T2]))
print(stats)
for e in ex:
    print("  ", e)

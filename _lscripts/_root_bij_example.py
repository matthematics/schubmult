import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from schubmult import RCGraph
from schubmult.combinatorics.inversions_tableau import InversionsTableau
p = 4
for rows in ([(4, 3, 2), (4,), (), ()], [(5, 3, 2), (4,), (), ()]):
    T = RCGraph(rows)
    print("T=", rows, "perm", T.perm, "hw?", all(T.raising_operator(i) is None for i in range(1, p)))
    comp = {(): T}; order = [()]; seen = {T: ()}; k = 0
    while k < len(order):
        path = order[k]; X = comp[path]; k += 1
        for i in range(1, p):
            Y = X.lowering_operator(i)
            if Y is not None and Y not in seen:
                seen[Y] = path + (i,); comp[path + (i,)] = Y; order.append(path + (i,))
    for path in order:
        X = comp[path]
        print("   f", path, [tuple(r) for r in X], sorted(InversionsTableau.from_rc_graph(X)._dict.items()))

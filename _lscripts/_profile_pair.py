import time
from schubmult import *
from schubmult.combinatorics.wc_graph import WCGraph
import _lscripts.squash_the_grove_wc_graphs as M

bw = M.bw
n = 5
length = n - 1
perm1 = uncode([0, 0, 0, 1])
perm2 = uncode([0, 0, 2, 1])


def T(label, fn):
    t = time.time()
    r = fn()
    print(f"{label}: {time.time() - t:.2f}s", flush=True)
    return r


t0 = time.time()
u1 = T("untagged perm1", lambda: M.untagged_groth_elem(perm1, length))
u2 = T("untagged perm2", lambda: M.untagged_groth_elem(perm2, length))
g1 = WCGraph.grove_wcs(perm1.pad_code(n - 1))
g2 = WCGraph.grove_wcs(perm2.pad_code(n - 1))


def poly_for(u, g, perm):
    return bw.from_dict({k: v for k, v in u.items() if bw.key_to_wc_graph(k).resize(n - 1) in g or bw.key_to_wc_graph(k).perm != perm})


p1 = T("poly_for perm1", lambda: poly_for(u1, g1, perm1))
p2 = T("poly_for perm2", lambda: poly_for(u2, g2, perm2))
print("len p1", len(p1), "len p2", len(p2), flush=True)
prod = T("poly1*poly2", lambda: p1 * p2)
print("len prod", len(prod), flush=True)
tp = T("to_wc_graph_ring_element", lambda: prod.to_wc_graph_ring_element())
print("len tp", len(tp), flush=True)
print(f"TOTAL: {time.time() - t0:.2f}s", flush=True)

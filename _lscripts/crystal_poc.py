from schubmult import *
from schubmult.rings.combinatorial.bounded_rc_factor_algebra import _elem_factor_from_rc
from schubmult.combinatorics.set_word import SetLetter

br = BoundedRCFactorAlgebra()
bw = BoundedWCFactorAlgebra()

def _set_from_lv(lv):
    """Return a set of integers from a length vector."""
    return {i for i, v in enumerate(lv, start=1) if v > 0}

def _set_to_lv(s, length):
    return tuple([1 if i in s else 0 for i in range(1, length + 1)])

def _setletter_to_elem_wc(setletter, p, k):
    if len(setletter) < p:
        return None
    return next(iter(WCGraph.all_wc_graphs(uncode([0]*(k - p) + [1] * p), k, weight=_set_to_lv(setletter, k))))
    

def test_crystal(n):
    perms = Permutation.all_permutations(n)
    for w in perms:
        the_wcs = set()
        for rc in RCGraph.all_rc_graphs(w, n - 1):
            rc_tensor = _elem_factor_from_rc(rc)
            base_p = [rc.perm.inv for rc in rc_tensor]
            setletter_tensor = CrystalGraphTensor(*[SetLetter(_set_from_lv(rc.length_vector), length=rc.perm.max_descent) for rc in rc_tensor])
            crystal = setletter_tensor.full_crystal
            for c in crystal:
                wc_graph_tensor_list = [_setletter_to_elem_wc(c.factors[i], base_p[i], rc_tensor[i].perm.max_descent) for i in range(len(c.factors))]
                if any(wc is None for wc in wc_graph_tensor_list):
                    continue
                wc_graph_tensor = bw.make_key(wc_graph_tensor_list, n - 1)
                wc_graph = bw.key_to_wc_graph(wc_graph_tensor)
                if wc_graph.perm != w:
                    continue#, f"Mismatch for {w}: {wc_graph.perm} != {w}\n{wc_graph=}\n{c=}\n{wc_graph_tensor=}\n{rc_tensor=}"
                the_wcs.add(wc_graph)
        assert the_wcs == set(WCGraph.all_wc_graphs(w, n - 1)), f"Mismatch for {w}, {len(the_wcs)} != {len(set(WCGraph.all_wc_graphs(w, n - 1)))}, {the_wcs=}\n{set(WCGraph.all_wc_graphs(w, n - 1))=}"
        print(f"Yay {w}")

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    test_crystal(n)
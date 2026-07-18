from schubmult import *
from schubmult.rings.polynomial_algebra import *

bw = BoundedWCFactorAlgebra()
br = BoundedRCFactorAlgebra()

def _wc_key_to_rc_key(key):
    return tuple([next(iter(RCGraph.elem_sym_rcs(len(wc.perm_word), wc.perm.max_descent, weight=wc.length_vector))) for wc in key])

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [perm.pad_code(n - 1) for perm in perms]
    for comp1, comp2 in itertools.product(comps, repeat=2):
        groth1 = bw.full_groth_elem(uncode(comp1), n, 1)
        groth2 = bw.full_groth_elem(uncode(comp2), n, 1)
        wc_to_factors = {bw.key_to_wc_graph(k).resize(n-1): k for k in set(groth1.keys()) | set(groth2.keys())}
        wc_to_factors = {wc: k for wc, k in wc_to_factors.items() if wc.perm.pad_code(n - 1) in {comp1, comp2}}
        prd = 0
        for wc1, wc2 in itertools.product(WCGraph.grove_wcs(comp1), WCGraph.grove_wcs(comp2)):#{rc for rc in RCGraph.all_rc_graphs(uncode(comp2), len(comp2)) if rc.forest_weight == comp2}):
            prd += bw(wc_to_factors[wc1]) * bw(wc_to_factors[wc2])
            #bw.from_dict({bw.make_key([WCGraph(k) for k in key], n): 1 for key, _ in br.from_rc_graph(br.key_to_rc_graph(br.make_key(_wc_key_to_rc_key(wc_to_factors[wc2]), n)), n).items()})
        prd = prd.to_wc_graph_ring_element().resize(n - 1)
        real_prd = GrovePoly(*comp1) * GrovePoly(*comp2)
        test_prd = 0
        for wc, coeff in prd.items():
            if wc.grove_weight == wc.length_vector and wc.grove_invariant == wc:
                test_prd += coeff * GrovePoly(*wc.grove_weight)
            #test_prd += coeff * GrovePoly(*wc.forest_weight)
        assert test_prd.almosteq(real_prd), f"Failed for {comp1}, {comp2}: {test_prd} != {real_prd}"
        print("Success for", comp1, comp2)
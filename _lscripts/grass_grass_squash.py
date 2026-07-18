from schubmult import *
from schubmult.rings.polynomial_algebra import *
from functools import cache

bw = BoundedWCFactorAlgebra()
wr = WCGraphRing()

def _wc_key_to_rc_key(key):
    return tuple([next(iter(RCGraph.elem_sym_rcs(len(wc.perm_word), wc.perm.max_descent, weight=wc.length_vector))) for wc in key])

@cache
def grass_poly(grass1, k):
    return wr.from_dict(dict.fromkeys(WCGraph.all_wc_graphs(grass1, k), 1))
    

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for k in range(1, n):
        permsk = [perm for perm in perms if perm.max_descent <= k]
        grass_permsk = [perm for perm in permsk if len(perm.descents()) == 1 and perm.max_descent == k]
        for normal1, grass2 in itertools.product(permsk, grass_permsk):
            grass_poly1 = grass_poly(normal1, k)
            grass_poly2 = grass_poly(grass2, k)
            result = wr.zero
            for wc1, coeff1 in grass_poly1.items():
                for wc2, coeff2 in grass_poly2.items():
                    result += coeff1 * coeff2 * wr(wc1.squash_product(wc2))

            #prd = (grass_poly1 * grass_poly2).to_wc_graph_ring_element().resize(k)
            
            real_prd = GrothendieckPoly(normal1, k) * GrothendieckPoly(grass2, k)
            test_prd = 0
            for wc, coeff in result.items():
                if wc.is_principal:
                    test_prd += coeff * GrothendieckPoly(wc.perm, k)
                #test_prd += coeff * GrovePoly(*wc.forest_weight)
            print("Testing", normal1, grass2)
            assert test_prd.almosteq(real_prd), f"Failed for {normal1}, {grass2}: {test_prd} != {real_prd}"
            print("Success for", normal1, grass2)
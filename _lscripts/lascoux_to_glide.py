from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import lascoux_poly

def _bwify_key(key):
    new_key = [next(iter(RCGraph.elem_sym_rcs(len(k.perm_word), k.perm.max_descent, weight=k.length_vector))) for k in key]
    return br.make_key(new_key, key.size)

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    for perm in perms:
        if perm.inv == 0:
            continue
        lascoux_by_hecke = {}
        the_min_by_hecke = {}
        
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            if wc.is_quasi_yamanouchi:# and wc.strong_hecke_invariant == principal.strong_hecke_invariant:
                lascoux_by_hecke[wc.strong_hecke_invariant] = lascoux_by_hecke.get(wc.strong_hecke_invariant, 0) + GlidePoly(*wc.length_vector)
                if wc.strong_hecke_invariant not in the_min_by_hecke or wc.length_vector < the_min_by_hecke[wc.strong_hecke_invariant]:
                    the_min_by_hecke[wc.strong_hecke_invariant] = wc.length_vector

        for hecke, poly in lascoux_by_hecke.items():
            min_key = the_min_by_hecke[hecke]
            actual_lascoux = GlidePoly.from_expr(lascoux_poly(min_key, Sx.genset), length=n-1)
            #result = poly.change_basis(GlidePolyBasis)
            assert poly.almosteq(actual_lascoux), f"Mismatch for {comp}: {poly} != {actual_lascoux}\n{poly-actual_lascoux}"
        #assert result.almosteq(actual_lascoux), f"Mismatch for {comp}: {result} != {actual_lascoux}\n{result-actual_lascoux}"
        print("Success for", perm)
            
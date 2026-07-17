from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import lascoux_poly
from schubmult.symbolic.poly.variables import *

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    full_snapped_mins = {}
    for perm in perms:
        for wc in WCGraph.all_wc_graphs(perm):
            if wc.is_quasi_yamanouchi:
                if wc.strong_hecke_invariant not in full_snapped_mins or wc.length_vector < full_snapped_mins[wc.strong_hecke_invariant]:
                    full_snapped_mins[wc.strong_hecke_invariant] = wc.length_vector

    for perm in perms:
        if perm.inv == 0:
            continue
        lascoux_by_hecke = {}
        the_min_by_hecke = {}
        
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            lascoux_by_hecke[wc.strong_hecke_invariant] = lascoux_by_hecke.get(wc.strong_hecke_invariant, set())
            lascoux_by_hecke[wc.strong_hecke_invariant] |= {wc}
            #lascoux_by_hecke.get(wc.strong_hecke_invariant, 0) + GlidePoly(*wc.length_vector)
            if wc.strong_hecke_invariant not in the_min_by_hecke or wc.length_vector < the_min_by_hecke[wc.strong_hecke_invariant]:
                the_min_by_hecke[wc.strong_hecke_invariant] = wc.length_vector

        for hecke, lascoux_set in lascoux_by_hecke.items():
            actual_lascoux = lascoux_poly(the_min_by_hecke[hecke], Sx.genset)
            for index in range(1, n - 1):
                hecke_poly = 0
                mgenset = MaskedGeneratingSet(Sx.genset, index_mask = tuple(range(1, index + 1)))
                MGlidePoly = PolynomialAlgebra(GlidePolyBasis(genset=mgenset))
                for wc in lascoux_set:
                    wc_upper, wc_lower = wc.vertical_cut(index)
                    wc_upper_normal = wc_upper.normalize()
                    wc_lower_normal = wc_lower.normalize()
                    #if wc_upper.is_quasi_yamanouchi and wc_lower.is_quasi_yamanouchi and 
                    if wc_upper_normal.length_vector == full_snapped_mins[wc_upper_normal.strong_hecke_invariant] and wc_lower_normal.length_vector == full_snapped_mins[wc_lower_normal.strong_hecke_invariant]:
                        #hecke_poly += PA.from_expr(lascoux_poly(wc_upper.length_vector, Sx.genset), length=n - 1) * PA.from_expr(lascoux_poly(wc_lower.length_vector, mgenset), length=n-1)
                        #hecke_poly += GlidePoly(*wc_upper.length_vector).expand() * MGlidePoly(*wc_lower.length_vector).expand()
                        hecke_poly += lascoux_poly(wc_upper.length_vector, Sx.genset) * lascoux_poly(wc_lower.length_vector, mgenset)
                #assert hecke_poly.almosteq(PA.from_expr(actual_lascoux, length=n-1)), f"Mismatch for {hecke}: {hecke_poly} != {actual_lascoux}\n{hecke_poly-actual_lascoux}"
                assert (hecke_poly - actual_lascoux).expand() == 0, f"Mismatch for {hecke}: {hecke_poly} != {actual_lascoux}\n{hecke_poly-actual_lascoux}"
        print("Success for", perm)
            
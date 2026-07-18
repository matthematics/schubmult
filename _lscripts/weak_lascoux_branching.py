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
    
    for perm in perms:
        if perm.inv == 0:
            continue
        lascoux_by_hecke = {}
        
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            lascoux_by_hecke[wc.strong_hecke_invariant] = lascoux_by_hecke.get(wc.strong_hecke_invariant, set())
            lascoux_by_hecke[wc.strong_hecke_invariant] |= {wc}
    
        for hecke, lascoux_set in lascoux_by_hecke.items():
            actual_lascoux = lascoux_poly(WCGraph._extremal_weight(perm, n - 1, hecke), Sx.genset)
            for index in range(1, n - 1):
                hecke_poly = 0
                mgenset = MaskedGeneratingSet(Sx.genset, index_mask = tuple(range(1, index + 1)))
                MGlidePoly = PolynomialAlgebra(GlidePolyBasis(genset=mgenset))
                for wc in lascoux_set:
                    wc_upper, wc_lower = wc.vertical_cut(index)
                    #if wc_upper.is_quasi_yamanouchi and wc_lower.is_quasi_yamanouchi and 
                    if wc_upper.length_vector == wc_upper.extremal_weight and wc_lower.length_vector == wc_lower.extremal_weight:
                        hecke_poly += lascoux_poly(wc_upper.length_vector, Sx.genset) * lascoux_poly(wc_lower.length_vector, mgenset)
                assert (hecke_poly - actual_lascoux).expand() == 0, f"Mismatch for {hecke}: {hecke_poly} != {actual_lascoux}\n{hecke_poly-actual_lascoux}"
        print("Success for", perm)
            
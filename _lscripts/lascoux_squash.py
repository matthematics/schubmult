from schubmult import *
from schubmult.rings.polynomial_algebra import *


if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.pad_code(n - 1)) for perm in perms]

    by_desc = {}
    for perm in perms:
        by_desc[perm.max_descent] = by_desc.get(perm.max_descent, set())
        by_desc[perm.max_descent] |= {perm}

    grass_by_desc = {}
    for desc, perms in by_desc.items():
        grass_by_desc[desc] = {perm for perm in perms if len(perm.descents()) == 1}

    wr = WCGraphRing()
    for k in range(1, n):
        for perm in by_desc[k]:
            comp = perm.pad_code(k)
            wcs = {wc for wc in WCGraph.all_wc_graphs(uncode(comp), k) if wc.strong_hecke_invariant == WCGraph.principal_wc(uncode(comp), k).strong_hecke_invariant}
            for grass_perm in grass_by_desc[k]:
                
                grass_poly = GrothendieckPoly(grass_perm, k).change_basis(LascouxPolyBasis)
                actual_result = LascouxPoly(*comp) * grass_poly
                try_result = 0
                
                for grass_wc in WCGraph.all_wc_graphs(grass_perm, k):
                    for wc in wcs:
                        squashed = wc.squash_product(grass_wc)
                        try_result += wr(squashed)

                invariant_weight = {}
                for wc, coeff in try_result.items():
                    if wc.strong_hecke_invariant not in invariant_weight:
                        invariant_weight[wc.strong_hecke_invariant] = wc.length_vector
                    elif invariant_weight[wc.strong_hecke_invariant] > wc.length_vector:
                        invariant_weight[wc.strong_hecke_invariant] = wc.length_vector

                result = 0
                for wc, coeff in try_result.items():
                    if wc.length_vector == invariant_weight[wc.strong_hecke_invariant]:
                        result += coeff * LascouxPoly(*wc.length_vector)
                assert result.almosteq(actual_result), f"Mismatch for {comp=} {grass_perm=} {k=}: {result} != {actual_result}\n{result-actual_result=}"
                print(f"Success for {comp=} {grass_perm=} {k=}: {result=}")
    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     lascoux_by_hecke = {}
    #     the_min_by_hecke = {}
        
    #     for wc in WCGraph.all_wc_graphs(perm, n - 1):
    #         if wc.is_quasi_yamanouchi:# and wc.strong_hecke_invariant == principal.strong_hecke_invariant:
    #             lascoux_by_hecke[wc.strong_hecke_invariant] = lascoux_by_hecke.get(wc.strong_hecke_invariant, 0) + GlidePoly(*wc.length_vector)
    #             if wc.strong_hecke_invariant not in the_min_by_hecke or wc.length_vector < the_min_by_hecke[wc.strong_hecke_invariant]:
    #                 the_min_by_hecke[wc.strong_hecke_invariant] = wc.length_vector

    #     for hecke, poly in lascoux_by_hecke.items():
    #         min_key = the_min_by_hecke[hecke]
    #         actual_lascoux = GlidePoly.from_expr(lascoux_poly(min_key, Sx.genset), length=n-1)
    #         #result = poly.change_basis(GlidePolyBasis)
    #         assert poly.almosteq(actual_lascoux), f"Mismatch for {comp}: {poly} != {actual_lascoux}\n{poly-actual_lascoux}"
    #     #assert result.almosteq(actual_lascoux), f"Mismatch for {comp}: {result} != {actual_lascoux}\n{result-actual_lascoux}"
    #     print("Success for", perm)
            
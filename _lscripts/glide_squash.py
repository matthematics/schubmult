from schubmult import *
from schubmult.rings.polynomial_algebra import *

def _bwify_key(key):
    new_key = [next(iter(RCGraph.elem_sym_rcs(len(k.perm_word), k.perm.max_descent, weight=k.length_vector))) for k in key]
    return br.make_key(new_key, key.size)

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    #comps = [tuple(perm.pad_code(n - 1)) for perm in perms]
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
            wcs = {wc for wc in WCGraph.all_wc_graphs(perm, k) if wc.dst == WCGraph.principal_wc(perm, k)}
            for grass_perm in grass_by_desc[k]:
                grass_poly = GlidePoly(*grass_perm.trimcode)
                actual_result = GlidePoly(*perm.trimcode) * grass_poly
                result = 0
                
                for grass_wc in WCGraph.all_wc_graphs(grass_perm, k):
                    if grass_wc.dst != WCGraph.principal_wc(grass_perm, k):
                        continue
                    for wc in wcs:
                        squashed = wc.squash_product(grass_wc)
                        #try_result += wr(squashed)
                        if squashed.is_quasi_yamanouchi:
                            result += GlidePoly(*squashed.perm.trimcode)

                assert result.almosteq(actual_result), f"Mismatch for {perm=} {grass_perm=} {k=}: {result} != {actual_result}\n{result-actual_result=}"
                print(f"Success for {perm=} {grass_perm=} {k=}: {result=}")
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
            
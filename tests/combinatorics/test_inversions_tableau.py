def test_groth_match():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.symbolic import S, expand, Symbol
    from schubmult.symbolic.poly.schub_poly import grothendieck_poly
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet
    from schubmult.combinatorics.inversions_tableau import InversionsTableau
    from schubmult.abc import x
    from schubmult import Permutation
    
    beta = Symbol("β")
    zz = ZeroGeneratingSet()

    for w in Permutation.all_permutations(4):
        result = S.Zero
        for wc_graph in WCGraph.all_wc_graphs(w, 3):
            result += InversionsTableau.from_wc_graph(wc_graph).polyvalue(x, beta=beta, prop_beta=True)
        via_it = result
        direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        assert expand(via_it - direct, deep=True) == S.Zero, f"Failed for {w}: {expand(via_it - direct, deep=True)=}"

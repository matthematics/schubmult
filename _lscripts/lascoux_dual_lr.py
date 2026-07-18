from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.symbolic.common_polys import lascoux_poly
from schubmult.symbolic.poly.variables import *

def _lascoux_graphs(comp, length):
    return {wc for wc in WCGraph.all_wc_graphs(uncode(comp), length) if wc.strong_hecke_invariant == WCGraph.principal_wc(uncode(comp), length).strong_hecke_invariant}

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    wr = WCGraphRing()
    LascouxDual = FreeAlgebra(LascouxBasis)
    for comp1, comp2 in itertools.product(comps, repeat=2):
        prd = wr(WCGraph.principal_wc(uncode(comp1), len(comp1))) * wr(WCGraph.principal_wc(uncode(comp2), len(comp2)))
        #print("Success for", comp1, comp2)
        real_prd = LascouxDual(*comp1) * LascouxDual(*comp2)
        test_prd = 0
        invariant_by_weight = {}
        for wc, coeff in prd.items():
            if wc.extremal_weight not in invariant_by_weight:
                invariant_by_weight[wc.extremal_weight] = wc.strong_hecke_invariant
            elif invariant_by_weight[wc.extremal_weight] != wc.strong_hecke_invariant:
                continue
            test_prd += coeff * LascouxDual(*wc.extremal_weight)
        assert test_prd.almosteq(real_prd), f"Failed for {comp1}, {comp2}: {test_prd} != {real_prd}"
        print(f"Success for {comp1}, {comp2}")
            
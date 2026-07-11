from schubmult import *  # noqa
from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.polynomial_algebra.forest_poly_basis import (
    ForestPolyBasis,
    _grove_polynomial_from_indfor,
)
from schubmult.symbolic import Symbol

beta = Gx._beta

Gx1 = GrothendieckRing(Gx.genset, beta=1)
ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))

r = RCGraphRing()
Grove1 = PolynomialAlgebra(GrovePolyBasis(Sx.genset, beta=1))

def grove_in_forests(comp):
    
    grothperms = Gx1.from_expr(Grove1(*comp).expand())
    # return ForestPoly.from_expr(grove)
    result = {}
    for permg, coeff in grothperms.items():
        
        for perm, coeff2 in WCGraph.groth_to_schub(permg, 1).items():
            forests = {rc for rc in RCGraph.all_rc_graphs(perm, len(comp)) if rc.is_forest_rc}
            for frc in forests:
                
                #if any(frc.length_vector in {wc.length_vector for wc in WCGraph.grove_wcs(comp)} for frc in forests):
                #print(f"{permg.pad_code(len(comp))=}, {coeff=}")
                result[frc.forest_weight] = result.get(frc.forest_weight, 0) + coeff * coeff2
            #elem_of_paint += , coeff*coeff2))
    
    # for stinkbat, coeff in elem_of_paint.items():
    # #for stinkbat in WCGraph.all_wc_graphs(uncode(comp), len(comp)):
    #     #if stinkbat.grove_weight == comp:
    #     if stinkbat.forest_weight == stinkbat.length_vector:
    #         result[stinkbat.forest_weight] = result.get(stinkbat.forest_weight, 0) + coeff
    return result

if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    perms = Permutation.all_permutations(n)
    comps = sorted({perm.pad_code(n - 1) for perm in perms})
    for comp in comps:
        if sum(comp) == 0:
            continue
        print(f"{comp=}")
        forest_real = ForestPoly.from_expr(Grove1(*comp).expand(), length=len(comp))
        expansion = ForestPoly.from_dict(grove_in_forests(comp))
        #forests = sorted(key for key, coeff in expansion.items() if coeff != 0)
        assert forest_real.almosteq(expansion), f"Mismatch for {comp}: {forest_real} != {expansion}"
        print(f"Good {expansion=}")
        #print(f"grove{comp} -> {forests}")

        

from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from sympy import pretty_print



if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    comps = [perm.pad_code(n - 1) for perm in Permutation.all_permutations(n)]
    
    ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    r = ForestRCGraphRing()
    #CompSchub = PolynomialAlgebra(CompositionSchubertBasis)
    length = n - 1
    for comp in comps:
        #for k in range(1, n):
        
        if True:
            k = length
            for p in range(1, k + 1):
                acc = r.zero
                elem_comp = tuple([0] * (k - p) + [1] * p)
                prodo = ForestPoly(*comp) * ForestPoly(*elem_comp)
                for forest_rc in RCGraph.all_forest_rcs(comp):
                    for elem_rc in RCGraph.all_rc_graphs(uncode(elem_comp), length):
                        the_squash = forest_rc.squash_product(elem_rc)
                        forest_squash = next(iter([rc for rc in RCGraph.all_forest_rcs(the_squash.forest_weight) if rc.length_vector == the_squash.length_vector and rc.omega_invariant[1] == the_squash.omega_invariant[1]]), None)
                        assert forest_squash.forest_weight in prodo
                        # acc[forest_squash.forest_weight] = acc.get(forest_squash.forest_weight, set())
                        # # if forest_squash in acc[forest_squash.forest_weight]:
                        # #     print(f"Duplicate squash product for {forest_rc} and {elem_rc} yielding {the_squash} with forest weight {the_squash.forest_weight}")
                        # acc[forest_squash.forest_weight].add(forest_squash)
                        acc += r(forest_squash)
                        #the_squash.forest_weight == the_squash.perm.pad_code(length), f"Squash product of {forest_rc} and {elem_rc} is {the_squash} and has forest weight {the_squash.forest_weight} instead of {the_squash.perm.pad_code(length)}"
                for frc, coeff in acc.items():
                    # if compo not in acc:
                    #     print(f"Missing squash product for {comp} and {elem_comp} yielding {compo}")
                    assert prodo[frc.forest_weight] == coeff, f"Coefficient mismatch for {comp} and {elem_comp} yielding {frc.forest_weight}, got {coeff} but expected {prodo[frc.forest_weight]}"
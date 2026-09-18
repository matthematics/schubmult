from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.rings.polynomial_algebra import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    r = RCGraphRing()
    Forest = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    comps = [tuple(perm.pad_code(n - 1)) for perm in Permutation.all_permutations(n)]
    for comp1, comp2 in itertools.product(comps, repeat=2):
        print(f"Trying {comp1}, {comp2}")
        forest_prod = Forest(*comp1) * Forest(*comp2)
        squash_result = 0
        for rc1, rc2 in itertools.product({rc for rc in RCGraph.all_rc_graphs(uncode(comp1), n - 1) if rc.forest_weight == comp1}, {rc for rc in RCGraph.all_rc_graphs(uncode(comp2), n - 1) if rc.forest_weight == comp2}):
            squash_result += r(rc1.squash_product(rc2))
        for rc, coeff in squash_result.items():
            if coeff != 0:
                assert forest_prod.get(rc.forest_weight) == coeff, f"Error: {rc.forest_weight} has coefficient {forest_prod.get(rc.forest_weight)} in the forest product but {coeff} in the squash product\n{forest_prod=}\n{squash_result=}"
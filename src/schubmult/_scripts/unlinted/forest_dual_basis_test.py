from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *

ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
ForestDual = FreeAlgebra(ForestBasis)

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    comps = [perm.pad_code(n-1) for perm in Permutation.all_permutations(n)]
    for comp1, comp2 in itertools.product(comps, repeat=2):
        Forest = ForestPoly(*comp1)
        dual = ForestDual(*comp2)
        if comp1 == comp2:
            assert dual.pairing(Forest) == 1, f"Failed for {comp1} and {comp2}, expected 1 but got {dual.pairing(Forest)}"
        else:
            assert dual.pairing(Forest) == 0, f"Failed for {comp1} and {comp2}, expected 0 but got {dual.pairing(Forest)}"
    print("Success")
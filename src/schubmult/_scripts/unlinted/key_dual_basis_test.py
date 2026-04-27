from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.rings.free_algebra import *

KeyPoly = PolynomialAlgebra(KeyPolyBasis(Sx.genset))
KeyDual = FreeAlgebra(KeyBasis)

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    comps = [perm.pad_code(n-1) for perm in Permutation.all_permutations(n)]
    for comp1, comp2 in itertools.product(comps, repeat=2):
        key = KeyPoly(*comp1)
        dual = KeyDual(*comp2)
        if comp1 == comp2:
            assert dual.pairing(key) == 1, f"Failed for {comp1} and {comp2}, expected 1 but got {dual.pairing(key)}"
        else:
            assert dual.pairing(key) == 0, f"Failed for {comp1} and {comp2}, expected 0 but got {dual.pairing(key)}"
    print("Success")
from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import lascoux_poly

if __name__ == "__main__":
    import sys
    import itertools
    from sympy import pretty_print
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        glides_by_grove = {}
        reps = set()
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            if wc.is_quasi_yamanouchi:
                glides_by_grove[wc.grove_invariant] = glides_by_grove.get(wc.grove_invariant, set())
                glides_by_grove[wc.grove_invariant] |= {wc}

        
        
        for grove_wc, reps in glides_by_grove.items():
            test_result = sum([GlidePoly(*wc.length_vector) for wc in reps]).change_basis(GrovePolyBasis)
            assert test_result.almosteq(GrovePoly(*grove_wc.length_vector)), f"Mismatch for {perm}: {test_result} != {GrovePoly(*grove_wc.length_vector)}\n{test_result-GrovePoly(*grove_wc.length_vector)}"
            print("Paint")
        
        print("Success for", perm)
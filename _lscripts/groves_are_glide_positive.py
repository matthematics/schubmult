from schubmult import *
from schubmult.rings.polynomial_algebra import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    #wr = WCGraphRing()
    for comp in comps:
        grove = GrovePoly(*comp)
        result = 0
        for wc in WCGraph.grove_wcs(comp):
            if wc.is_quasi_yamanouchi:
                result += GlidePoly(*wc.length_vector)
        assert result.change_basis(GrovePolyBasis).almosteq(grove), f"Mismatch for {comp}: {result.change_basis(GrovePolyBasis)} != {grove}"
        print("Success for", comp)
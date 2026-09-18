from schubmult import *
from schubmult.rings.polynomial_algebra import PolynomialAlgebra, ForestPolyBasis, Schub, SchubertPolyBasis


if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    non_forest_patterns = [[1,4,3,2],[2,4,1,3],[2,4,3,1],[1,4,5,2,3],[3,2,1,5,4],[3,4,1,2,6,5]]
    forest_good = []
    r = BoundedRCFactorAlgebra()
    for perm in perms:
        if all(not perm.has_pattern(pat) for pat in non_forest_patterns):
            #print(perm.trimcode)
            forest_good.append(perm)
    # ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    #for perm1, perm2 in itertools.product(forest_good, forest_good):
    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm2.has_pattern([1,4,3,2]) or perm2.has_pattern([3,1,2]):
            continue
        elem1 = r.schub_elem(perm1, n)
        elem2 = r.schub_elem(perm2, n)
        result = 0
        for key1, coeff1 in elem1.items():
            rc1 = r.key_to_rc_graph(key1)
            for key2, coeff2 in elem2.items():
                rc2 = r.key_to_rc_graph(key2)
                prd = r.from_rc_graph(rc1,n) * r.from_rc_graph(rc2,n)
                result += coeff1 * coeff2 * prd.to_rc_graph_ring_element()
    
        pants = 0
        expected = Sx(perm1) * Sx(perm2)
        for rc, coeff in result.items():
            if rc.is_principal or (rc.perm not in forest_good and rc.forest_weight == rc.length_vector):
                pants += coeff * Sx(rc.perm)
        if expected != pants:
            print(f"Error: Mismatch for {perm1.trimcode} * {perm2.trimcode}: expected {expected}, got {pants}")
            exit(1)
        else:
            print(f"Success: {perm1.trimcode} * {perm2.trimcode} = {pants}")
from schubmult import *
from schubmult.rings.polynomial_algebra import *


def pieri_forest(n):
    """???"""
    from schubmult.rings.combinatorial.forest_rc_ring import _canonical_rc
    perms = Permutation.all_permutations(n)
    # seen = {}
    
    r = RCGraphRing()
    for k in range(n, n + 1):
        for p in range(1, k + 1):
            for perm in perms:
                pieris = {}
                for rc in RCGraph.all_rc_graphs(perm, n):
                    for elem_sym_rc in RCGraph.all_rc_graphs(uncode([0] * (k - p) + [1] * p), n):
                        pieris[(rc.forest_invariant,p)] = pieris.get((rc.forest_invariant, p), 0) + r(rc.squash_product(elem_sym_rc))
                    #stump_key
                for (forest_inv, degree), rc_elem in pieris.items():
                    for rcc, coeff in rc_elem.items():
                        assert all(coeff == v for rc2, v in rc_elem.items() if rcc.forest_invariant == rc2.forest_invariant), f"Forest invariant {forest_inv} has multiple RC graphs with different coefficients: {rc_elem}"
                        assert all([(rc in rc_elem.keys()) for rc in RCGraph.all_rc_graphs(rcc.perm, n) if rc.forest_invariant == rcc.forest_invariant]), f"Forest invariant {forest_inv} is missing RC graphs: {rc_elem}"
                    print("Good forest invariant", forest_inv, "with degree", degree)
                        # print(f"Forest invariant {forest_inv} has

    


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    # cauchy_dual_key(n)
    # cauchy_dual_forest(n)
    pieri_forest(n)